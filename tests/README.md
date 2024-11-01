# Tests
- Unit test: Checks basic BMI functionality and runs model for one timestep (1 hour) to compare results against a benchmark.
- Synthetic tests for LGAR: There are three synthetic examples for testing and to demonstrate the usage of the model. The example simulates 12 hours of rainfall, infiltration, soil saturation, surface ponding, and surface runoff. The examples use different soil hydraulic properties and different precipitation intensities without PET to simulate precipitation partitioning into infiltration and runoff. These simulate precipitation partitionind due to a short rainfall pulse that infiltrates entirely, followed by a longer precipitation pulse that generates runoff.
- Synthetic tests for LGARTO: in the process of merging LGARTO and LGAR, we have added 4 synthetic tests that demonstrate LGARTO's ability to emulate Richards equation solver (HYDRUS) results. These are 72 hour long tests with a varying number of precipitation pulses (2 or 0), AET pulses (3 or 0), over a variety of layered soil scenarios, intended to demonstrate the functionality of LGARTO in a detailed way. 
- Tests using real forcing and soils data from USDA SCAN sites for LGARTO: also in the process of merging LGARTO to the main branch, we have added 4 year long simulations that use real soils and forcing data that show one one could practically expect LGARTO to compare against Richards equation results. Currently both synthetic and real LGARTO simulations are run from the main LGAR-C folder (not tests).
- Note that, in LGARTO's development, we have identified opportunities to speed up the code while reducing its accuracy compared to the Richards equation results. Accordingly, most of the LGARTO runs here have been generated using slight modifications to LGAR.cxx, where: the number of initial wetting fronts was increased from 4 to 16, wetting fronts are inserted near the bottom of the root zone when the difference between capillary head of adjacent wetting fronts is smaller, and new wetting fronts are more frequently inserted at the soil surface. Essentially, inserting TO wetting fronts more frequently such that there are more wetting fronts that are closer in capillary head values can modestly increase accuracy at a substantial runtime cose. Finally, note that setting use_closed_form_G to true in the config file will improve model speed but decrease the accuracy of the capillary drive component in both the calculation for dzdt for surface wetting fronts, as well as initial wetting front depth. 
- Comparisons of LGARTO and LGAR: to ensure that the development work resulting in LGARTO did not affect LGAR, a notebook has been added demonstrating the lack of significant changes between LGARTO with TO mode set to off, and LGAR, instances of 4 simulations (synth1, synth2, Bushland, and Phillipsburg).
- Stability testing: Because LGARTO and LGAR aim to emulate the Rchards equation, whose solution can be rather complex, a large number of edge cases that cause LGAR or LGARTO to crash (via infinite loops, large mass balance errors, or other misc. problems) are possible. Throughout the development of these models, we have constantly tested boradly for stability, running the models with a variety of randomly generated parameter sets over multiple forcing datasets. We have added notebooks that determine the frequency of unstable model runs. While both LGARTO and LGAR currently have 100% stability with the provided stability tests, there is always a chance that new edge cases could present themselves. 
- Animation code: often the best way to investigate a LGAR / LGARTO model run is to visualize how the wetting fronts move over time, together with the effect on cumulative fluxes. Code that provides such visualizations has been added. 

## Unittest
### Build
```
mkdir build && cd build (inside LGAR-C directory)
cmake ../ -DUNITTEST=ON
make && cd ../tests
```

### Run:
run `./run_unittest.sh`

## Synthetic tests
### Build
```
mkdir build && cd build (inside LGAR-C directory; if build already exists then clean it)
cmake ../ -DSTANDALONE=ON
make && cd ../tests
```

### Run:
run `./run_synthetic.sh OPTION` (for synthetic test; OPTION = 1 or 2 - these numbers correspond to different examples)


#### Visualization
  - Use `plot_synthetic_examples.ipynb` to plot and compare synthetic lgar examples with hydrus output

#### The unittest performs several checks:
  1. Check names, number, and memory allocation of BMI input/output variables
  2. Check units of BMI input/output variables
  3. Check number of layers and number of wetting fronts against the benchmark
  4. Check BMI Grid function (e.g., grid_id, grid_size, etc.)
  5. Test `GetVar*` methods for the BMI input variables and compare against initial (or prescribed) data; if failed, will throw an error
  6. Loop over the input variables, use `Set*` and `Get*` methods to verify `Get*` return the same data set by `Set*`
  7. Using `Update` method, advance the model to get updated depths and soil moisture of the wetting fronts. Compare against the benchmark values.

  #### Unit test results
  If everything goes well, you should see the following

  $\textcolor{green}{\text{| ************************************************************} }$ \
  $\textcolor{green}{\text{| All BMI Tests passed: YES} }$ \
  $\textcolor{green}{\text{| LASAM Calibration test = YES} }$ \
  $\textcolor{green}{\text{| ************************************************************} }$

