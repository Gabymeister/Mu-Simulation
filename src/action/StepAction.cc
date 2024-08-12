/*
Tom: this step action is used to only record K0/KS/KL decays

*/

#include "action.hh"
#include <tls.hh>
#include <G4MTRunManager.hh>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
namespace MATHUSLA { namespace MU {


TTree* StepAction::_step_data = nullptr;
std::string StepAction::file_name = ".root";
std::string StepAction::prefix = ".";
bool StepAction::step_data_valid = false;
// TFile * StepAction::f = new TFile(StepAction::file_name.c_str(), "UPDATE");

double StepDataStore::_step_x = 0;
double StepDataStore::_step_y = 0;
double StepDataStore::_step_z = 0;
double StepDataStore::_step_x_end = 0;
double StepDataStore::_step_y_end = 0;
double StepDataStore::_step_z_end = 0;
double StepDataStore::_step_px = 0;
double StepDataStore::_step_py = 0;
double StepDataStore::_step_pz = 0;
double StepDataStore::_energy_loss = 0;
double StepDataStore::_deposit = 0;
int StepDataStore::_pdg = 0;
int StepDataStore::_trackid = 0;
int StepDataStore::_trackid_parent = 0;
int StepDataStore::_trackid_status = 0;
int StepDataStore::_material_index = 0;
int StepDataStore::_evt_number = 0;


void StepDataStore::Initialize(int id){
    _evt_number = id;

    auto event_id = std::to_string(id);
    auto tree_name = "event_" + event_id;
    StepAction::step_data_valid = false;

    StepAction::_step_data = new TTree(tree_name.c_str(),tree_name.c_str());

    StepAction::_step_data->Branch("X_S", &_step_x, "X_S/D");
    StepAction::_step_data->Branch("Y_S", &_step_y, "Y_S/D");
    StepAction::_step_data->Branch("Z_S", &_step_z, "Z_S/D");
    StepAction::_step_data->Branch("X_END_S", &_step_x_end, "X_END_S/D");
    StepAction::_step_data->Branch("Y_END_S", &_step_y_end, "Y_END_S/D");
    StepAction::_step_data->Branch("Z_END_S", &_step_z_end, "Z_END_S/D");    

    StepAction::_step_data->Branch("PX_S", &_step_px, "PX_S/D");
    StepAction::_step_data->Branch("PY_S", &_step_py, "PY_S/D");
    StepAction::_step_data->Branch("PZ_S", &_step_pz, "PZ_S/D");

    StepAction::_step_data->Branch("DE", &_energy_loss, "DE/D");
    StepAction::_step_data->Branch("DEPOSIT", &_deposit, "DEPOSIT/D");
    StepAction::_step_data->Branch("PDG", &_pdg, "PDG/I");
    StepAction::_step_data->Branch("TRACK", &_trackid, "TRACK/I");
    StepAction::_step_data->Branch("TRACK_PARENT", &_trackid_parent, "TRACK_PARENT/I");
    StepAction::_step_data->Branch("TRACK_STATUS", &_trackid_status, "TRACK_STATUS/I");
    StepAction::_step_data->Branch("MATERIAL", &_material_index, "MATERIAL/I");
    StepAction::_step_data->Branch("EVENT_NUMBER", &_evt_number, "EVENT_NUMBER/I");
  }

//__Action Initialization Constructor___________________________________________________________
StepAction::StepAction() : G4UserSteppingAction() {
  // std::cout<<"-------------------------------------------------"<<std::endl;
  // std::cout<<file_name<<std::endl;
  // std::cout<<"-------------------------------------------------"<<std::endl;
  // f =  TFile::Open(file_name.c_str(), "RECREATE");
}
//----------------------------------------------------------------------------------------------

void StepAction::UserSteppingAction( const G4Step* step){
  

  const auto step_point = step->GetPreStepPoint();
  const auto post_step_point = step->GetPostStepPoint();
  const auto position   = G4LorentzVector(step_point->GetGlobalTime(), step_point->GetPosition());
  const auto position_end   = G4LorentzVector(post_step_point->GetGlobalTime(), post_step_point->GetPosition());
  const auto momentum   = G4LorentzVector(step_point->GetTotalEnergy(), step_point->GetMomentum());

  auto pdg = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();

  // Check if it is KL/KS, or the secondaries contains KL/KS
  std::set<int> KL_pdgIDs {130, 310, 311, -311};
  bool is_KLKS=false;
  bool is_KLKS_secondary=false;

  if(KL_pdgIDs.count(pdg) != 0){
    is_KLKS = true;
  }
  auto secondary = step->GetSecondaryInCurrentStep();
  size_t size_secondary = secondary->size();
  if (size_secondary){
    for (size_t i=0; i<(size_secondary);i++){
      auto secstep = (*secondary)[i];
      int particle_pdg = secstep->GetParticleDefinition()->GetPDGEncoding();
      if(KL_pdgIDs.count(particle_pdg) != 0)
        is_KLKS_secondary=true;
    }
  }   
  if (!(is_KLKS || is_KLKS_secondary))
    return;

  StepAction::step_data_valid = true;

  StepDataStore::_step_x = position.x();
  StepDataStore::_step_y = position.y();
  StepDataStore::_step_z = position.z();
  StepDataStore::_step_x_end = position_end.x();
  StepDataStore::_step_y_end = position_end.y();
  StepDataStore::_step_z_end = position_end.z();  

  StepDataStore::_step_px = momentum.x();
  StepDataStore::_step_py = momentum.y();
  StepDataStore::_step_pz = momentum.z();

  StepDataStore::_deposit = step->GetTotalEnergyDeposit();

  StepDataStore::_energy_loss = step_point->GetKineticEnergy() - post_step_point->GetKineticEnergy();

  StepDataStore::_pdg = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
  StepDataStore::_trackid = step->GetTrack()->GetTrackID();
  StepDataStore::_trackid_parent = step->GetTrack()->GetParentID();
  StepDataStore::_trackid_status = step->GetTrack()->GetTrackStatus();

  StepDataStore::_material_index = step_point->GetMaterial()->GetIndex();
  _step_data->Fill();

}

void StepAction::WriteTree(int id){
  auto event_id = std::to_string(id);
  // std::cout<<"Write step to file: "<<file_name<<std::endl;
  if (step_data_valid ){
    TFile* f = new TFile(file_name.c_str(), "UPDATE");
    f->cd();
    _step_data->Write();
    f->Write();
    f->Close();

    delete f;
    delete _step_data;

    // TFile f(file_name.c_str(), "UPDATE");
    // f.cd();
    // _step_data->Write();
    // f.Write();
    // f.Close();    
  }

}

} } /* namespace MATHUSLA::MU */


