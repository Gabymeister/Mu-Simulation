# MATHUSLA Mu-Simulation


_simulation of muons through earth material using Geant4_


## Build & Run

The simulation comes with a simple build script called `install` which allows for build customization and execution of the muon simulation.

Here is a list of useful commands:

| Action             | Options for `./install` |
|:------------------:|:-----------------------:|
| Build Only         | `(none)`                |
| Build with CMake   | `--cmake` (use this for the first time)|
| Build and Auto Run | `--run`                 |
| Clean CMake Build  | `--cmake --clean`       |
| More Options       | `--help`                |

After building, the executable, `simulation`, is moved into the root directory of the project.

The simulation executable itself comes with several configuration parameters:

| Action                | Short Options    | Long Options        |
|:---------------------:|:----------------:|:-------------------:|
| Event Count           | `-e <count>`     | `--events=<count>`  |
| Particle Generator    | `-g <generator>` | `--gen=<generator>` |
| Detector              | `-d <detector>`  | `--det=<detector>`  |
| Custom Script         | `-s <file>`      | `--script=<file>`   |
| Data Output Directory | `-o <dir>`       | `--out=<dir>`       |
| Number of Threads     | `-j <count>`     | `--threads=<count>` |
| Visualization         | `-v`             | `--vis`             |
| Save All Generator Events         | `NA` | `--save_all`        |
| Save Events With Pseudo-Digi Cuts (only for Cosmic geometry) | `NA` | `--cut_save`        |
| Bias Muon Nuclear interaction in Earth (for Cosmic and Box geometry) | `NA` | `--bias`        |
| Turn On Five Body Muon Decays     | `-f` | `--five_muon`       |
| Non-Random Five Body Decays       | `-n` | `--non_random`      |
| Debug Mode (write out step files)  | `-D`             | `--debug`            |
| Quiet Mode            | `-q`             | `--quiet`           |
| Help                  | `-h`             | `--help`            |

Note: The Five Body Muon Decays option will only save tracks with a five-body decay in a certain zone in the detector. Be sure this is what you want.

**1. Running the simulation interactively**

The simulation exectuable can be run interactively using the -v option. For example, the following command will run 1 event with the Box detector.

```
./simulation -v -e 1
```

This is not recommanded unless for debugging purposes, becuase there are many setting commands that needs to be entered later by hand.

**2. Running the simulation with script**

The simulation exectuable can take a script that list all the setting and commands to be executed.

```
./simulation -d DETECTOR_TYPE  -s YOUR_SCIPT.mac
```

Available commands are summarized in later section. Here is an example of a script:

```
# example1.mac

# Select the Detector. Box is the full 100x100m tracker
# Although selected here, you still need to add '-d Box' option when running the simulation
/det/select Box

# Select the Generator
/gen/select basic
# Set Generator-specific particle properties
/gen/basic/id 13
/gen/basic/t0 0 ns
/gen/basic/vertex 0   0   84.57 m
/gen/basic/p_unit 0 0 -1
/gen/basic/p_mag 100 GeV/c
/gen/basic/phi 0 rad
/gen/basic/eta  1

# Set number of evets to run
/run/beamOn 100
```

Arguments can also be passed through the simulation to a script. Adding key value pairs which correspond to aliased arguments in a script, will be forwarded through. Here's an example to make the momentum (p_mag) adjustable with arguments. Replace the number 100 in the line for momentum with a key name {momentum}:  `/gen/range/p_mag {momentum} GeV/c` where `{ }` denotes a key. The script will be evaluated by substituting the key for its value from the argument. Run the simulation with the following command:

```
./simulation -s example1.mac momentum 100
```

This will pass the key value pairs `(momentum 100)` to the underlying script `example1.mac` and set momentum to 100.

A general command to run the simulation is like this:

```
./simulation -q  -o OUTPUT_DIR  -d DETECTOR_TYPE -s YOUR_SCRIPT.mac key1 100 key2 10 ....
```

**IMPORTANT**
Even if you set the detector type in your script, it is still necessary to specify it using '-d DETECTOR_TYPE' option when running the simulation.

### Generators

| Generator  | Details                                               |
|:----------:|:-----------------------------------------------------:|
| basic         | Generate particles with FIXED position and momentum |
| range         | Generate particles with position and momentum in a range|
| polar         | Generate particles with position and momentum in a range, in polar coordinates|
| pythia        | similar to the range generator|
| file_reader   | Generate particles with user-specified position and momentum for **EACH** event, can be used as a vertex gun|

There are three general purpose generators built in, `basic`, `range`, and `polar`. The `basic` generator produces a particle with constant `pT`, `eta`, and `phi` while the `range` generator produces particle within a specified range of values for each of the three variables. Any variable can also be fixed to a constant value. The `polar` generator uses the angles from spherical coordinates, polar and azimuth, along with an energy input to generate particles. The polar angle in `polar` generator can be either a constant or within a specified range, while the azimuth is only given within a range.

There is also a _Pythia8_ generator installed which behaves similiarly to the `range` generator.

There is a `file_reader` generator that produces particles with properties that are specified in an input file. The `file_reader` does not do randomization. Each entry in the input file correspond to exactly one event and may include multiple particles, which means that the physics process to generate the primary partiles needs to be taken care of outside GEANT4. This generator provide the possibility to generate controllable primary vertex at given location with pre-assigned momentum for each secondary particles.

The generator defaults are specified in `src/action/GeneratorAction.cc` but they can be overwritten by a custom generation script.

Examples of generators can be found in `scripts/generators/`

### Custom Detector

A custom Detector can be specified at run time from one of the following installed detectors:

| Detector   | Status    | Details                                               |
|:----------:|:---------:|:-----------------------------------------------------:|
| Prototype  | COMPLETED | Test stand for the MATHUSLA project                   |
| Box        | COMPLETED | Large Detector as seen in MATHUSLA Original Schematic |
| Cosmic     | COMPLETED | Identical detector design as in Box but optimized for cosmic studies |
| Flat       | BUILDING  | Cheaper Alternative to Box                            |
| MuonMapper | COMPLETED | Measures Muon Energies after Rock Propagation         |

### Custom Scripts

A custom _Geant4_ script can be specified at run time. The script can contain generator specific commands and settings as well as _Pythia8_ settings in the form of `readString`. The script can also specify the detector to use during the simulation. This section summarize all available commands.

**1. Selecting detector**

    /det/select Box,Prototype,Cosmic,Flat,MuonMapper

**2. Selecting Generator**

    /gen/select basic,range,polar,pythia,file_reader

**3. Generator-specific settings**

Basic

| Command    | Usage    |
|:----------:|:---------:|
| /gen/basic/id  | COMPLETED |

**4. General Geant4 commands**

/run/beamOn 100

## Digitizer
Included in this repository is a program which operates on the output of the main simulation. It combines truth hits that fall within an indistinguishable spacetime window of each other (the parameters of which are specified by the user). It also adds dark counts, the rate of which may also be specified by the user. Further details are described in the digitizer directory.

### Installation
The digitizer may be installed with the equivalent script `install_digitizer`. When installing for the first time, be sure to use the option `--cmake`.

## VectorExtraction
This directory contains scripts relevant to the generation and formatting of MG5 simulations, which are the primary source for studying collider-generated muons. 
### MadGraph and Pythia8 ### 

This version of the code uses MadGraph v3.5.1. MadGraph is a self-contained framework for simulating SM and BSM phenomenology, computing cross-sections, and generating hard-events. Pythia8 may also be installed in MadGraph, which simulates multi-state parton showers.

Reference:

    J. Alwall et al, "The automated computation of tree-level and next-to-leading order differential cross sections, and their matching to parton shower simulations"

    Bierlich, Christian, et al. "A comprehensive guide to the physics and usage of PYTHIA 8.3." SciPost Physics Codebases (2022): 008.

Usage: 
Madgraph may be downloaded from this link: http://launchpad.net/madgraph5/3.0/3.5.x/+download/MG5_aMC_v3.5.4.tar.gz

Unzip it:

    gzip -d MG5_aMC_v3.5.4.tar.gz

Then you can run it using the command

    ./MG5_aMC_v3_5_4/bin/mg5_aMC

This should provide an interactive session. The first time this is done, run the command

    install pythia8

Then type "exit." This will permanently install pythia8 into the self-contained MadGraph framework.

It is recommended before running the program to set the appropriate paths, as certain cross-sections may conflict between your default pythia version and MadGraph's version. One should set the following path before running MadGraph with Pythia8:

    export PYTHIA8=/project/def-mdiamond/tomren/mathusla/pythia8308
    export PYTHIA8DATA="MG5_DIR/HEPTools/pythia8/share/Pythia8/xmldoc"

where `MG5_DIR` is the MadGraph Directory (so from before, `MG5_DIR=MG5_aMC_v3_5_4`)

### Running MadGraph ### 
One may execute a script for MadGraph using the command 

    MG5_dir/bin/mg5_aMC your_script.txt

Examples of these scripts are found in `VectorExtraction/MadGraphScripts`. 

The output of MadGraph is not able to be read immediately by the file_reader generator, and must be formatted. The script `run_muon_extract.py` operates on the hepmc file generated by a madgraph run and outputs a text file in the format appropriate for the file_generator. It requires three arguments: The first is the hepmc file. The second is the desired location of the output. The third is the minimum pT of events that are desired. 

`combine_muon_data.py` takes several output files from `run_muon_extract.py` and combines them into a single file while also generating the macro for Geant4 to run. It takes the following four inputs: The location where the final output will be stored, the ID for this simulation (i.e. the job ID), the number of files to be combined, and the path to the files that are to be combined. 

It is recommended to use the scripts in the Mu-Background repository which submit jobs that run both these operations. 

The output is a ROOT file. In the file, there are some metadata and a ROOT TTree. The tree name depends on the detector follwing the convention of \'{detectortype}_run\'. For example, a box detector will have tree name "box_run".

Quantities in the ROOT TTree:

| Key | Data type | Unit |Comment|
|:-------|:----|:----:|:----|
|Hit_x,y,z|vector| \[cm\]| In CMS coordinate|
