## Configuration File
Example configuration files are provided in this directory. To build and run the given examples see the instructions [here](https://github.com/NOAA-OWP/LGAR-C/blob/master/INSTALL.md#build-standalone-mode).

A detailed description of the parameters for model configuration (i.e., initialize/setup) is provided below. 


| Variable | Datatype |  Limits  | Units | Role | Process | Description |
| -------- | -------- | ------ | ----- | ---- | ------- | ----------- |
| forcing_file | string | - | - | filename | - | provides precip. and PET inputs |
| soil_params_file | string | - | - | filename | - | provides soil types with van Genuchton parameters |
| layer_thickness | double (1D array)| - | cm | state variable | - | individual layer thickness (not absolute)|
| initial_psi | double (scalar)| >=0 | cm | capillary head | - | used to initialize layers with a constant head |
| ponded_depth_max | double (scalar)| >=0 | cm | maximum surface ponding | - | the maximum amount of water unavailable for surface drainage, default is set to zero |
| timestep | double (scalar)| >0 | sec/min/hr | temporal resolution | - | timestep of the model |
| forcing_resolution | double (scalar)| - | sec/min/hr | temporal resolution | - | timestep of the forcing data. Recommended value of 3600 seconds. |
| endtime | double (scalar)| >0 | sec, min, hr, d | simulation duration | - | time at which model simulation ends |
| layer_soil_type | int (1D array) | - | - | state variable | - | layer soil type (read from the database file soil_params_file) |
| max_soil_types | int | >1 | - | - | - | maximum number of soil types read from the file soil_params_file (default is set to 15) |
| wilting_point_psi | double (scalar) | - | cm | state variable | - | wilting point (the amount of water not available for plants) used in computing AET. Suggested value is 15495.0 cm, corresponding to 15 atm. |
| field_capacity_psi | double (scalar) | - | cm | state variable | - | capillary head corresponding to volumetric water content at which gravity drainage becomes slower, used in computing AET. Suggested value is 340.9 cm for most soils, corresponding to 1/3 atm, and 103.3 cm for sands, corresponding to 1/10 atm. |
| use_closed_form_G | bool | true or false | - | - | - | determines whether the numeric integral or closed form for G is used; a value of true will use the closed form. This defaults to false. |
| giuh_ordinates | double (1D array)| - | - | state parameter | - | GIUH ordinates (for giuh based surface runoff) |
| verbosity | string | high, low, none | - | debugging | - | controls IO (screen outputs and writing to disk) |
| sft_coupled | Boolean | true, false | - | model coupling | impacts hydraulic conductivity | couples LASAM to SFT. Coupling to SFT reduces hydraulic conducitivity, and hence infiltration, when soil is frozen|
| soil_z | double (1D array) | - | cm | spatial resolution | - | vertical resolution of the soil column (computational domain of the SFT model) |
| calib_params | Boolean | true, false | - | calibratable params flag | impacts soil properties | If set to true, soil `smcmax`, `smcmin`, `vg_n`, `vg_alpha`, `hydraulic_conductivity`, `field_capacity_psi`, and `ponded_depth_max` are calibrated. defualt is false. vg = van Genuchten, SMC= soil moisture content |
| adaptive_timestep | Boolean | true, false | - | adaptive timestep flag | impacts timestep | If set to true, LGAR will use an internal adaptive timestep, and the above timestep is used as a minimum timestep (recommended value of 300 seconds). The adaptive timestep will never be larger than the forcing resolution. If set to false, LGAR will use the above specified timestep as a fixed timestep. Testing indicates that setting this value to true substantially decreases runtime while negligibly changing the simulation. We recommend this to be set to true. |
