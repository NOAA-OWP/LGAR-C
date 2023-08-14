### Model setup/run

#### Build
 - To build LASAM (see [instructions](https://github.com/NOAA-OWP/LGAR-C/blob/master/INSTALL.md))

#### Run (in LGAR-C directory)
```
mkdir lasam && cd lasam
ln -s ../extern
ln -s ../data
cp extern/LGAR-C/data/vG_default_params.dat data/
../cmake_build/ngen data/catchment_data.geojson cat-27 data/nexus_data.geojson nex-26 extern/LGAR-C/realizations/realization_config_lasam.json
```