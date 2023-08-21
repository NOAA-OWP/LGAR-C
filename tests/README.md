# Tests
- Unit test: Checks basic BMI functionality and runs model for one timestep (1 hour) to compare results against a benchmark.
- Synthetic tests: There are three synthetic examples for testing and to demonstrate the usage of the model. The example simulates 12 hours of rainfall, infiltration, soil saturation, surface ponding, and surface runoff. The examples use different soil hydraulic properties and different precipitation intensities without PET to simulate precipitation partitioning into infiltration and runoff. These simulate precipitation partitionind due to a short rainfall pulse that infiltrates entirely, followed by a longer precipitation pulse that generates runoff.

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
run `./run_synthetic.sh OPTION` (for synthetic test; OPTION = 0, 1, or 2 - these numbers correspond to different examples)


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

