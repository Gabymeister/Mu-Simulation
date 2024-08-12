#include <iostream>
#include "RunManager.hh"
#include "TreeHandler.hh"
#include "NoiseMaker.hh"
#include "Digitizer.hh"
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "par_handler.hh"
#include <map>

int RunManager::StartTracking()
{
// ---------- Declaring Tree Handlers and setting up geometry ---------------------------------------
	TH = new TreeHandler(_InputTree_Name, _InputFile_Name, _OutputTree_Name, OutFileName());
	if (TH->IsNull()) {
		std::cout << "Sorry, I couldn't open that file" << std::endl;
		return 0;
	} 
	
	TreeHandler *CH;
	if (cosmic) {
		CH = new TreeHandler(_InputTree_Name, _CosmicFile_Name, TH);
		if (CH->IsNull()) {
			std::cout << "Sorry, I couldn't open that cosmic file" << std::endl;
			return 0;
		}
	} 

	GeometryHandler _geometry(_geomTree_Name, _InputFile_Name);
	_geometry.LoadGeometry();

	int events_handled = 0;

	int dropped_hits = 0;
	int floor_wall_hits = 0;

	//Total LHC muons
	int nLHC = 0;
	//Number of total noise events
	int nNoise = 0;
	//Number of cosmic events
	int totalCosmic = 0;
	//Full Total
	int FullNumber = 0;

	ParHandler hndlr;
	hndlr.Handle();

// ---------- Finding excluded layers (if any) ------------------------------------------------------
	//Getting the layers that we will not be including
	std::vector<int> excludedLayers;
	int diff = (int)(_geometry.GetYLayers().size()/2) - hndlr.par_map["n_layers"];
	for (int i = 0; i < diff; i++){
		// Remove the horizontal tracking layers
		int layerNum = (int)(_geometry.GetZWalls().size()/2) 
		+ (int)(_geometry.GetYFloor().size()/2) 
		+ i;
		excludedLayers.push_back(layerNum);
		// Remove the back wall tracking layers
		layerNum = (int)(_geometry.GetZWalls().size()/2) 
		+ (int)(_geometry.GetYFloor().size()/2) 
		+ (int)(_geometry.GetYLayers().size()/2) 
		+ (int)(_geometry.GetZBack().size()/2)
		- i - 1;
		excludedLayers.push_back(layerNum);
	}
		
// ---------- Loop through main tree ----------------------------------------------------------------
	if (hndlr.par_map["branch"] == 1.0) std::cout << "Running in Cosmic Mode" << std::endl;

	_digitizer->par_handler = &hndlr;

	NoiseMaker::preDigitizer(&_geometry);
	
	_digitizer->InitGenerators();
	
	while (TH->Next() >= 0)
	{
		if (events_handled >= hndlr.par_map["end_ev"])//cuts::end_ev)
			break;
		if (events_handled >= hndlr.par_map["start_ev"]) //cuts::start_ev)
		{
			if (hndlr.par_map["debug"] == 1) 
				std::cout << "\n"<< std::endl;
			if ((events_handled) % 1000 == 0 || hndlr.par_map["debug"] == 1 || hndlr.par_map["debug_vertex"] == 1 )
				std::cout << "=== Event is " << events_handled <<" ==="<< std::endl;

			TotalEventsProcessed++;
			_digitizer->clear();
			TH->LoadEvent();

			int n_added = 0; //This tracks the number of sim_hits added, needed for setting proper index
			int n_TH = 0; // Number of sim_hits added from TH
			//adding all hits of the tree into the digitizer
			for (int n_hit = 0; n_hit < TH->sim_numhits; n_hit++){
				//physics::sim_hit *current = new physics::sim_hit(TH, &_geometry, n_hit);
				physics::sim_hit *current = new physics::sim_hit(TH, &_geometry, n_hit);
				//Check if sim hit is in an excluded layer
				current->index = n_added;
				TH->sim_hit_type_buf->push_back(0);
				current->SetType(0);
				if (hndlr.par_map["branch"] == 1.0) {
					current->x += detector::COSMIC_SHIFT[0];
					current->y += detector::COSMIC_SHIFT[1];
					current->z += detector::COSMIC_SHIFT[2];
				}
				if (excludedLayers.size() == 0 || std::find(excludedLayers.begin(), excludedLayers.end(), current->det_id.layerID) 
				== excludedLayers.end()) { //Not an exluded layer, so we add it
					_digitizer->AddHit(current);
					n_TH++;
				}
				nLHC++;
				n_added++;
			}

			//Doing cosmic (if any)
			//The following is used to assign unique track IDs that would otherwise be duplicated
			int n_cosmic = 0;	
			std::map<int, int> trackIDChanger;
			int maxID;
			int curID;
			if (cosmic){
				n_cosmic = _digitizer->generator.Poisson(hndlr.par_map["cosmic_rate"]);
				maxID = 0; // For setting new track IDs so there are no duplicates
				for (int n_hit = 0; n_hit < n_TH; n_hit++){
					if (_digitizer->hits[n_hit]->track_id > maxID) {
						maxID = _digitizer->hits[n_hit]->track_id;
					}
				}
				maxID++;
			}
			//adding cosmic events as chosen by poisson distribution
			while (n_cosmic > 0) {
				CH->index =(int)_digitizer->generator.Uniform(CH->NumEntries);
				if (CH->Next() < 0) {
					CH->index = 0;// for now, just loop back to start	
				}
				CH->LoadEvent();
				for (int n_hit = 0; n_hit < CH->sim_numhits; n_hit++){
					physics::sim_hit *current = new physics::sim_hit(CH, &_geometry, n_hit);
					//if (excludedLayers.size() != 0 && std::find(excludedLayers.begin(), excludedLayers.end(), current->det_id.layerID) 
					//!= excludedLayers.end()) { 
						//std::cout << "Layer ID " << current->det_id.layerID << std::endl;
				//		continue;
				//	}
					current->index = n_added;
					CH->sim_hit_type_buf->push_back(1);
					current->SetType(1);
					//Reassigning track IDs to new ones
					curID = current->track_id;
					for (int oldHit = 0; oldHit < n_TH; oldHit++) {
						if (curID == _digitizer->hits[oldHit]->track_id) {
							//This one not added yet
							if (trackIDChanger.find(curID) == trackIDChanger.end()){
								trackIDChanger[curID] = maxID;
								current->track_id = maxID;
								maxID++;
							} else {
								current->track_id = maxID;
							}
						}
					}
					if (hndlr.par_map["branch"] == 1.0) {
						current->x += detector::COSMIC_SHIFT[0];
						current->y += detector::COSMIC_SHIFT[1];
						current->z += detector::COSMIC_SHIFT[2];
					}
					if (excludedLayers.size() == 0 || std::find(excludedLayers.begin(), excludedLayers.end(), current->det_id.layerID) 
					== excludedLayers.end()) { // Then it's not excluded so we digitize it
						_digitizer->AddHit(current);
					}
					n_added++;
				}
				totalCosmic++;
				n_cosmic -= 1;
			}
			std::vector<physics::digi_hit *> digi_list = _digitizer->Digitize();
			std::vector<physics::digi_hit*> noise_digis;
			if(NoiseMaker::run){
				NoiseMaker* noise = new NoiseMaker(digi_list);
				noise_digis = noise->return_digis();
				for(auto digi:noise_digis){
					digi->SetType(2);
					digi_list.push_back(digi);
					nNoise++;
                		}
        		}
			TH->ExportDigis(digi_list, _digitizer->seed);
			TH->Fill();

			dropped_hits += _digitizer->dropped_hits;
			floor_wall_hits += _digitizer->floor_wall_hits;
		}
		events_handled++;
	}
	FullNumber = nLHC + nNoise + totalCosmic;
	std::cout << "Number of LHC events: " << nLHC << std::endl;
	std::cout << "Number of noise events: " << nNoise <<std::endl;
	std::cout << "Number of cosmic events: " << totalCosmic << std::endl;
	std::cout << "Number of events total: " << FullNumber << std::endl;
	if (hndlr.file_opened) {
		std::cout << "writing to file" << std::endl;
		TH->Write();
	}
	return 0;
}
