# Build and Run Instructions

Detailed instructions on how to build and run LASAM in two modes (standalone and nextgen framework). Building LASAM requires [GCC](https://gcc.gnu.org) and [CMAKE](https://cmake.org/) on your machine.

## Standalone mode example
The examples provided here simulate infiltration and surface runoff using the data from field sites located in Phillipsburg, Kansas, and Bushland, Texas. 
### Build
 - mkdir build && cd build (inside LGAR-C directory)
 - cmake ../ -DSTANDALONE=ON
 - make && cd ..
### Run
```
./build/lasam_standalone configs/config_lasam_X.txt (X = Phillipsburg, Bushland; run from LGAR-C directory)
```

## Nextgen framework example
See general [instructions](https://github.com/NOAA-OWP/ngen/wiki/NGen-Tutorial#running-cfe) for building models in the nextgen framework. Assuming you have a running nextgen framework, follow the below instructions to build LASAM and SLoTH, and then run the example.
### Build
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

### Run
```
 mkdir lasam && cd lasam (inside nextgen directory)
 ln -s ../extern
 ln -s ../data
 cp extern/LGAR-C/data/vG_default_params.dat data/
 ../cmake_build/ngen data/catchment_data.geojson cat-27 data/nexus_data.geojson nex-26 extern/LGAR-C/realizations/realization_config_lasam.json
```

### LASAM Coupling to Soil Freeze Thaw (SFT) Model
- Follow the instructions on [SoilFreezeThaw](https://github.com/NOAA-OWP/SoilFreezeThaw) repo to build SFT model and soil moisture profiles. Note [SoilMoistureProfiles](https://github.com/NOAA-OWP/SoilMoistureProfiles) is needed for the coupling of LASAM to SFT.
- Realization file for LASAM-SFT coupling is provided in the [realizations](./realizations/) directory (realization_config_lasam_sft.json)

  **Note:** Make sure the `"library_file"` and `"init_config"` in the BMI blocks in the realization files are pointing to the right files.
