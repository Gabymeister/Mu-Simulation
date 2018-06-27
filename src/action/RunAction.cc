/*
 * src/action/RunAction.cc
 *
 * Copyright 2018 Brandon Gomes
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "action.hh"

#include <fstream>
#include <ostream>

#include <Geant4/G4Threading.hh>
#include <Geant4/G4AutoLock.hh>
#include <Geant4/tools/wroot/file>

#include <ROOT/TFile.h>
#include <ROOT/TNamed.h>

#include "analysis.hh"
#include "geometry/Construction.hh"
#include "physics/Units.hh"

#include "util/io.hh"
#include "util/time.hh"
#include "util/stream.hh"

namespace MATHUSLA { namespace MU {

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Static Variables____________________________________________________________________________
G4ThreadLocal std::string _data_dir{};
G4ThreadLocal std::string _path{};
G4ThreadLocal uint_fast64_t _event_count{};
G4ThreadLocal uint_fast64_t _run_count{};
static G4Mutex _mutex = G4MUTEX_INITIALIZER;
//----------------------------------------------------------------------------------------------

//__Write Entry to ROOT File____________________________________________________________________
template<class... Args>
void _write_entry(const std::string& name,
                  Args&& ...args) {
  std::stringstream stream;
  util::stream::forward(stream, args...);
  TNamed entry(name.c_str(), stream.str().c_str());
  entry.Write();
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__RunAction Constructor_______________________________________________________________________
RunAction::RunAction(const std::string& data_dir) : G4UserRunAction() {
  _data_dir = data_dir == "" ? "data" : data_dir;
}
//----------------------------------------------------------------------------------------------

//__Run Initialization__________________________________________________________________________
void RunAction::BeginOfRunAction(const G4Run* run) {
  _path = _data_dir;
  util::io::create_directory(_path);
  _path += '/' + util::time::GetDate();
  util::io::create_directory(_path);
  _path += '/' + util::time::GetTime();
  util::io::create_directory(_path);
  _path += "/run" + std::to_string(_run_count) + ".root";

  _event_count = run->GetNumberOfEventToBeProcessed();

  Analysis::ROOT::Setup();
  Analysis::ROOT::Open(_path);
  Analysis::ROOT::CreateNTuple(
    Construction::Builder::GetDetectorDataName(),
    Construction::Builder::GetDetectorDataKeys());
}
//----------------------------------------------------------------------------------------------

//__Post-Run Processing_________________________________________________________________________
void RunAction::EndOfRunAction(const G4Run*) {
  if (!_event_count)
    return;

  Analysis::ROOT::Save();

  G4AutoLock lock(&_mutex);
  if (!G4Threading::IsWorkerThread()) {
    TFile root_file(_path.c_str(), "UPDATE");
    if (!root_file.IsZombie()) {
      root_file.cd();
      _write_entry("FILETYPE", "MATHULSA MU-SIM DATAFILE");
      _write_entry("DET", Construction::Builder::GetDetectorName());
      for (const auto& entry : GeneratorAction::GetGenerator()->GetSpecification())
        _write_entry(entry.name, entry.text);
      _write_entry("RUN", _run_count);
      _write_entry("EVENTS", _event_count);
      _write_entry("TIMESTAMP", util::time::GetString("%c %Z"));
      root_file.Close();
    }
    ++_run_count;
    std::cout << "\nEnd of Run\nData File: " << _path << "\n\n";
  }
  lock.unlock();
}
//----------------------------------------------------------------------------------------------

//__Get Current RunID___________________________________________________________________________
size_t RunAction::RunID() {
  return _run_count;
}
//----------------------------------------------------------------------------------------------

//__Get Current EventCount______________________________________________________________________
size_t RunAction::EventCount() {
  return _event_count;
}
//----------------------------------------------------------------------------------------------

} } /* namespace MATHUSLA::MU */
