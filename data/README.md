## Description of calibratable parameters

A detailed description of the calibratable parameters is given below. 

| Parameter name | Units |  Physical limits  | Range tested for stability | Applies to individual soil layers or entire model domain | Description |
| --- | --- | --- | --- | --- | --------------- |
| theta_r | - | 0< theta_r <1, <br> theta_r < theta_e | 0.01 < theta_r < 0.15 | Soil layer | theta_r is the residual water content, or the minimum volumetric water content that a soil layer can naturally attain. Note that theta_r must be less than than theta_e. This is set per soil layer, in the .dat file in the data directory.|
| theta_e | - | 0 < theta_e < 1, <br> theta_r < theta_e | 0.3 < theta_e < 0.8 | Soil layer | theta_e is the maximum volumetric water content that a soil layer can naturally attain. Note that theta_e must be greater than theta_r. This is set per soil layer, in the .dat file in the data directory.|
| alpha | 1/cm | alpha > 0 | 0.001 < alpha < 0.3 | Soil layer | alpha is the van Genuchten parameter related to the inverse of air entry pressure. Note that in nature, expected values of alpha are distributed logarithmically, so calibrating on the log of alpha rather than alpha directly is likely a better choice for most calibration algorithms. This is set per soil layer, in the .dat file in the data directory.|
| n | - | n > 1 | 1.01 < n <3 | Soil layer | n is the van Genuchten parameter related to pore size distribution. Note that another commonly reported van Genuchten parameter, m, is related to n via m = 1 - 1/n. We recommend that n > 1.01. Values larger than 3 are acceptable. This is set per soil layer, in the .dat file in the data directory.|
| Ks | cm/h | Ks > 0 | 0.001 < K_s < 100 | Soil layer | Ks is the saturated hydraulic conductivity of a soil. Note that in nature, expected values of Ks are distributed logarithmically, so calibrating on the log of Ks rather than Ks directly is likely a better choice for most calibration algorithms. This is set per soil layer, in the .dat file in the data directory.|
| ponded_head_max | cm | ponded_head_max >= 0 | 0 <= ponded_head_max <= 5 | Entire model domain | This is the maximum amount of ponded water that is allowed to accumulate on the soil surface. While stability tests have only included a maximum value of 5 cm, any value greater than or equal to 0 should be acceptable. A common choice will be 0. This parameter can be set in the config file.  |
| field_capacity_psi | cm | 0 < field_capacity_psi, <br>field_capacity_psi < wilting_point_psi | 10.3 < field_capacity_psi < 516.6 | Entire model domain | This is the wilting point of the model domain, expressed as a capillary head. Together with wilting_point_psi, the field capacity is used to determine the intensity of the reduction of PET to become AET. The numbers 10.3 cm and 516.6 cm correspond to pressures of 1/100 atm and 1/2 atm of water. Note that the model generally uses absolute values of capillary head; in this case, these limits are absolute values of negative numbers and physically represent unsaturated soil. While field capacity will vary per soil type, we use a single value for the entire model domain, following the method for PET->AET correction used by HYDRUS. This parameter can be set in the config file. |

Parameters that are specified per soil layer can either be scalar values with double precision (in the event the model is run with 1 layer) or vectors of doubles (in the event that the model is run with more than 1 layer), whereas parameters that are specified for the entire model domain are scalars with double precision. 

Stability testing efforts have included varying initial conditions, where the initial condition throughout the soil moisture profile is set with a single value for capillary head, ranging from 3000 cm to 10 cm (again note that these are absolute values of negative numbers and indicate unsaturated soils). This value can be set in the config file. 

We recommend that practical parameter values for calibration efforts of real soils use ranges that are for some parameters somewhat more restricted than the ones used in stability testing. For example, while it was desirable to test theta_e values up to 0.8, we do not expect to often see theta_e values of this magnitude in nature. Further, randomly sampling parameter values within even restricted ranges to build parameter sets could theoretically yield unrealistic soils (for example, a parameter set could have a K_s value of 0.01 cm/h, indicative of clay, and an n value of 2.5, indicative of sand). 

Below is a table of parameters for soils from the HYDRUS soils catalog, which can give insights to likely parameter values per soil class. It contains example parameters for 12 soil classes, taken from the HYDRUS-1D soils catalog, which in turn are based on the paper: Carsel, R.F., and Parrish, R. S., Developing joint probability distributions of soil water retention characteristics, Water Resour. Res, 24, 755-769, 1988. θr, θs, and n are unitless, α has units of 1/cm, and Ks has units of cm/h.

| Textural Class      | θr    | θs  | α     | n    | Ks          |
|---------------------|-------|-----|-------|------|-------------|
| Sand                | 0.045 | 0.43 | 0.145 | 2.68 | 29.7        |
| Loamy Sand          | 0.057 | 0.41 | 0.124 | 2.28 | 14.59166667 |
| Sandy Loam          | 0.065 | 0.41 | 0.075 | 1.89 | 4.420833333 |
| Loam                | 0.078 | 0.43 | 0.036 | 1.56 | 1.04        |
| Silt                | 0.034 | 0.46 | 0.016 | 1.37 | 0.25        |
| Silty Loam          | 0.067 | 0.45 | 0.02  | 1.41 | 0.45        |
| Sandy Clay Loam     | 0.1   | 0.39 | 0.059 | 1.48 | 1.31        |
| Clay Loam           | 0.095 | 0.41 | 0.019 | 1.31 | 0.26        |
| Silty Clay Loam     | 0.089 | 0.43 | 0.01  | 1.23 | 0.07        |
| Sandy Clay          | 0.1   | 0.38 | 0.027 | 1.23 | 0.12        |
| Silty Clay          | 0.07  | 0.36 | 0.005 | 1.09 | 0.02        |
| Clay                | 0.068 | 0.38 | 0.008 | 1.09 | 0.2         |

