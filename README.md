## Lumped Arid/Semi-arid Model (LASAM) for infiltration and surface runoff
The LASAM simulates infiltration and runoff based on Layered Green & Ampt with redistribution (LGAR) model. LGAR is a model which partitions precipitation into infiltration and runoff, and is designed for use in arid or semi arid climates. LGAR closely mimics precipitation partitioning results simulated by the famous Richards/Richardson equation (RRE), without the inherent reliability and stability challenges the RRE poses. Therefore, this model is useful when accurate, stable precipitation partitioning simulations are desired in arid or semi arid areas. LGAR has its python version too that is available [here](https://github.com/NOAA-OWP/LGAR-Py).

**Published papers:** For details about the model please see our manuscript on LGAR [link](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022WR033742)

### Instructions for: 
  - Installation (see [instructions](https://github.com/NOAA-OWP/LGAR-C/blob/master/INSTALL.md))
  - Test examples
    - simulations with synthetic forcing data, and unittest (see [tests](https://github.com/NOAA-OWP/LGAR-C/blob/master/tests/README.md))
    - simulations with real forcing data (see [tests](https://github.com/NOAA-OWP/LGAR-C/blob/master/INSTALL.md))
  - LASAM coupling to Soil Freeze Thaw (SFT) model (see [instructions](https://github.com/NOAA-OWP/LGAR-C/blob/master/INSTALL.md))

### Model Configuration File
  - Detailed description of the parameters for model configuration is provided ([here](https://github.com/NOAA-OWP/LGAR-C/tree/master/configs/README.md))
  
### Getting help
For questions, please contact Ahmad (ahmad.jan@noaa.gov) and/or Peter (peter.lafollette@noaa.gov), the two main developers/maintainers of the repository.

### Known issues or raise an issue
LASAM is a newly developed model and we are constantly looking to improve the model and/or fix bugs as they arise. Please see the Git Issues for known issues or if you want to suggest adding a capability or found a bug, please open an issue.
