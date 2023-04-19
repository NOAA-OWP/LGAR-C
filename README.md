## Lumped Arid/Semi-arid Model (LASAM) for infiltration and surface runoff
The LASAM simulates infiltration and runoff based on Layered Green & Ampt with redistribution (LGAR) model. LGAR is a model which partitions precipitation into infiltration and runoff, and is designed for use in arid or semi arid climates. LGAR closely mimics precipitation partitioning results simulated by the famous Richards/Richardson equation (RRE), without the inherent reliability and stability challenges the RRE poses. Therefore, this model is useful when accurate, stable precipitation partitioning simulations are desired in arid or semi arid areas. LGAR has its python version too that is available [here](https://github.com/NOAA-OWP/LGAR-Py).

**Published papers:** (provide a link here once LGAR paper is published)

### Build and Run LASAM
Here are two examples to build LASAM: 1) standalone mode and 2) ngen framework. Building LASAM requires [GCC](https://gcc.gnu.org) and [CMAKE](https://cmake.org/) on your machine.

**Example description:** The example provided here simulates infiltration and surface runoff for the period `2016-10-01 to 2017-08-09` using the data from a field site located in Phillipsburg, Kansas. There is also another example using data from Bushland, Texas. 

**Additional examples:**
- LASAM unit test. It is recommended to build and run the unittest before running other examples. More instructions are provided [here](https://github.com/NOAA-OWP/LGAR-C/tree/master/tests).
- Synthetic simulation ([here](https://github.com/NOAA-OWP/LGAR-C/tree/master/tests)) for testing/educational purposes. It simulates 12 hours of rainfall, infiltration, soil saturation, surface ponding, and surface runoff.


### LASAM standalone example
#### Build
 - git clone https://github.com/NOAA-OWP/LGAR-C
 - cd LGAR-C && mkdir build && cd build
 - cmake ../ -DSTANDALONE=ON
 - make && cd ..
#### Run
```
./build/lasam_standalone configs/config_lasam_Phillipsburg.txt (run this from LGAR-C directory)
```
### LASAM ngen framework example
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

## Configuration File
Example configuration files for the two examples above are provided in the [configs](./configs/) directory contains. The parameters in those configuration files are explained in the Table below.

| Variable | Datatype |  Limits  | Units | Role | Process | Description |
| -------- | -------- | ------ | ----- | ---- | ------- | ----------- |
| forcing_file | string | - | - | filename | - | provides precip. and PET inputs |
| soil_params_file | string | - | - | filename | - | provides soil types with van Genuchton parameters |
| layer_thickness | double (1D array)| - | cm | state variable | - | individual layer thickness (not absolute)|
| initial_psi | double (scalar)| >=0 | cm | capillary head | - | used to initialize layers with a constant head |
| ponded_depth_max | double (scalar)| >=0 | cm | maximum surface ponding | - | the maximum amount of water unavailable for surface drainage, default is set to zero |
| timestep | double (scalar)| >0 | sec/min/hr | temporal resolution | - | timestep of the model |
| forcing_resolution | double (scalar)| - | sec/min/hr | temporal resolution | - | timestep of the forcing data |
| endtime | double (scalar)| >0 | sec, min, hr, d | simulation duration | - | time at which model simulation ends |
| layer_soil_type | int (1D array) | - | - | state variable | - | layer soil type (read from the database file soil_params_file) |
| max_soil_types | int | >1 | - | - | - | maximum number of soil types read from the file soil_params_file (default is set to 15) |
| wilting_point_psi | double (scalar) | - | cm | state variable | - | wilting point (the amount of water not available for plants) used in computing AET |
| use_closed_form_G | bool | true or false | - | - | - | determines whether the numeric integral or closed form for G is used; a value of true will use the closed form. This defaults to false. |
| giuh_ordinates | double (1D array)| - | - | state parameter | - | GIUH ordinates (for giuh based surface runoff) |
| verbosity | string | high, low, none | - | debugging | - | controls IO (screen outputs and writing to disk) |
| sft_coupled | Boolean | true, false | - | model coupling | impacts hydraulic conductivity | couples LASAM to SFT. Coupling to SFT reduces hydraulic conducitivity, and hence infiltration, when soil is frozen|
| soil_z | double (1D array) | - | cm | spatial resolution | - | vertical resolution of the soil column (computational domain of the SFT model) |
