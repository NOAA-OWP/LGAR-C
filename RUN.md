# Running instructions
Here, we provide a few examples to run and test LASAM.

### Examples (standalone mode)
  - Unit test: It is recommended to build and run the unittest before running other examples. More instructions are provided [here](https://github.com/NOAA-OWP/LGAR-C/tree/ajk/doc_update/tests).
  - Synthetic simulation: ([here](https://github.com/NOAA-OWP/LGAR-C/tree/ajk/doc_update/tests)) for testing/educational purposes. It simulates 12 hours of rainfall, infiltration, soil saturation, surface ponding, and surface runoff.
 - Real fied example: The examples provided here simulates infiltration and surface runoff using the data from field sites located in Phillipsburg, Kansas and Bushland, Texas. 

### Run
 - Unittest ([here](https://github.com/NOAA-OWP/LGAR-C/tree/ajk/doc_update/tests))
 - Synthetic test ([here](https://github.com/NOAA-OWP/LGAR-C/tree/ajk/doc_update/tests))
 - Run real field examples (standalone model)
 ```
 ./build/lasam_standalone configs/config_lasam_Phillipsburg.txt (run from LGAR-C directory)
 ```
 - Run real field examples (nextgen mode)
 ```
 mkdir lasam && cd lasam (in LGAR-C directory)
 ln -s ../extern
 ln -s ../data
 cp extern/LGAR-C/data/vG_default_params.dat data/
 ../cmake_build/ngen data/catchment_data.geojson cat-27 data/nexus_data.geojson nex-26 extern/LGAR-C/realizations/realization_config_lasam.json
 ```
