# Installation instructions

Detailed instructions on how to install, configure, and run LASAM simulation examples.

### Build and Run LASAM
Here are two examples to build LASAM: 1) standalone mode and 2) ngen framework. Building LASAM requires [GCC](https://gcc.gnu.org) and [CMAKE](https://cmake.org/) on your machine.

**Example description:** The example provided here simulates infiltration and surface runoff for the period `2016-10-01 to 2017-08-09` using the data from a field site located in Phillipsburg, Kansas. There is also another example using data from Bushland, Texas. 

**Additional examples:**
- LASAM unit test. It is recommended to build and run the unittest before running other examples. More instructions are provided [here](https://github.com/NOAA-OWP/LGAR-C/tree/master/tests).
- Synthetic simulation ([here](https://github.com/NOAA-OWP/LGAR-C/tree/master/tests)) for testing/educational purposes. It simulates 12 hours of rainfall, infiltration, soil saturation, surface ponding, and surface runoff.

### LASAM Standalone Example
#### Build
 - git clone https://github.com/NOAA-OWP/LGAR-C
 - cd LGAR-C && mkdir build && cd build
 - cmake ../ -DSTANDALONE=ON
 - make && cd ..
 
#### Run
```
./build/lasam_standalone configs/config_lasam_Phillipsburg.txt (run from LGAR-C directory)
```

### LASAM Nextgen Framework Example
- See general [instructions](https://github.com/NOAA-OWP/ngen/wiki/NGen-Tutorial#running-cfe) for building models in the ngen framework. Assuming you have a running ngen framework, follow the below instructions to build LASAM and SLoTH, and then run the example.
#### Build

- #### LASAM
   - cd extern
   - git clone https://github.com/NOAA-OWP/LGAR-C (this should be removed when LGAR-C becomes a subrepo of the framework)
   - cmake -B extern/LGAR-C/cmake_build -S extern/LGAR-C/ -DNGEN=ON
   - make -C extern/LGAR-C/cmake_build/

- #### SLoTH
   SLoTH is also needed to run LASAM in the ngen framework. SLoTH is a BMI that is used to set a bmi variable(s) that is not provided by other BMIs but required by the model. So build [SLoTH](https://github.com/NOAA-OWP/SLoTH) using the following instructions
   - cd extern/sloth/ && git checkout latest
   - git submodule update --init --recursive
   - cd ../..
   - cmake -B extern/sloth/cmake_build -S extern/sloth/
   - make -C extern/sloth/cmake_build

#### Run
```
mkdir lasam && cd lasam
ln -s ../extern
ln -s ../data
cp extern/LGAR-C/data/vG_default_params.dat data/
../cmake_build/ngen data/catchment_data.geojson cat-27 data/nexus_data.geojson nex-26 extern/LGAR-C/configs/realization_config_lasam.json
```

### LASAM Coupling to Soil Freeze Thaw (SFT) Model
- Follow the instructions on [SoilFreezeThaw](https://github.com/NOAA-OWP/SoilFreezeThaw) repo to build SFT model and soil moisture profiles. Note [SoilMoistureProfiles](https://github.com/NOAA-OWP/SoilMoistureProfiles) is needed for the coupling of LASAM to SFT.
- Realization file for LASAM-SFT coupling is provided in the [configs](./configs/) directory (realization_config_lasam_sft.json)