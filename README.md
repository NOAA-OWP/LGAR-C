# Lumped Arid/Semi-arid Model (LASAM) for infiltration and surface runoff
The LASAM simulates infiltration and runoff based on Layered Green & Ampt with redistribution (LGAR) model. LGAR is a model which partitions precipitation into infiltration and runoff, and is designed for use in arid or semi-arid climates. LGAR closely mimics precipitation partitioning results simulated by the famous Richards/Richardson equation (RRE), without the inherent reliability and stability challenges the RRE poses. Therefore, this model is useful when accurate, stable precipitation partitioning simulations are desired in arid or semi-arid areas. LGAR in Python (no longer supported) is available [here](https://github.com/NOAA-OWP/LGAR-Py).

LASAM is designed for use in environments where cumulative potential evapotranspiration is greater than cumulative precipitation. Because the lower boundary condition of LASAM is effectively no-flow, the model assumes that water only leaves the vadose zone via AET. This is a reasonable assumption in arid and semi arid areas. If applied in humid areas, the model domain of LASAM will likely become completely saturated.

**Published papers:** For details about the model please see our manuscript on LGAR ([weblink](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022WR033742)).

## Build and Run Instructions
Detailed instructions on how to build and run LASAM can be found here [INSTALL](https://github.com/NOAA-OWP/LGAR-C/blob/master/INSTALL.md).
- Test examples highlights
  - simulations with synthetic forcing data and unittest (see [build/run](https://github.com/NOAA-OWP/LGAR-C/blob/master/tests/README.md)). 
  - simulations with real forcing data (see [build/run](https://github.com/NOAA-OWP/LGAR-C/blob/master/INSTALL.md#standalone-mode-example))
  - LASAM coupling to Soil Freeze Thaw (SFT) model (see [instructions](https://github.com/NOAA-OWP/LGAR-C/blob/master/INSTALL.md#lasam-coupling-to-soil-freeze-thaw-sft-model))

## Model Configuration File
A detailed description of the parameters for model configuration is provided [here](https://github.com/NOAA-OWP/LGAR-C/tree/master/configs/README.md).

## Calibratable parameters
A detailed description of calibratable parameters is provided [here](https://github.com/NOAA-OWP/LGAR-C/tree/master/data/README.md).

## Nextgen Realization Files
Realization files for running LASAM (coupled/uncoupled modes) in the nextgen framework are provided [here](https://github.com/NOAA-OWP/LGAR-C/tree/master/realizations/README.md).
  
## Getting help
For questions, please contact Ahmad (ahmad.jan(at)noaa.gov) and/or Peter (peter.lafollette(at)noaa.gov), the two main developers/maintainers of the repository.

## Known issues or raise an issue
LASAM is a newly developed model and we are constantly looking to improve the model and/or fix bugs as they arise. Please see the Git Issues for known issues or if you want to suggest adding a capability or to report a bug, please open an issue.

## Getting involved
See general instructions to contribute to the model development ([instructions](https://github.com/NOAA-OWP/LGAR-C/blob/master/CONTRIBUTING.md)) or simply fork the repository and submit a pull request.
