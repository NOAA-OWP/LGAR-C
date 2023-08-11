## Configuration File
A detailed description of the parameters for model configuration (i.e., initialize/setup) is provided below. The example configuration files are provided in the [configs](./configs/) directory.


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
| calib_params | Boolean | true, false | - | calibratable params flag | impacts soil properties | If set to true, soil `smcmax`, `vg_m`, and `vg_alpha` are calibrated. defualt is false |