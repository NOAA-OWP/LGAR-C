## Tests
#### - Unit test: Checks basic BMI functionality and runs model for one timestep (1 hour) to compare results against a benchmark.
#### - Synthetic test: For testing/educational purposes. It simulates 12 hours of rainfall, infiltration, soil saturation, surface ponding, and surface runoff.

#### Build:
  - mkdir build && cd build (inside LGAR-C directory)
  - cmake ../ -DUNITTEST:BOOL=ON (for unittest)
  - cmake ../ -DSTANDALONE:BOOL=ON (for synthetic test)
  - make

#### Run:
  - cd test
  - run `./run_unittest.sh` (for unittest)
  - run `./run_synthetic.sh` (for synthetic test)

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
  $\textcolor{green}{\text{| All BMI Tests passed: Yes} }$ \
  $\textcolor{green}{\text{| Infiltration_mm: (benchmark vs computed) | 1.896 vs 1.896} }$ \
  $\textcolor{green}{\text{| PET_mm: (benchmark vs computed) | 0.104 vs 0.104} }$ \
  $\textcolor{green}{\text{| AET_mm: (benchmark vs computed) | 0.0109182 vs 0.0109182} }$ \
  $\textcolor{green}{\text{| ************************************************************} }$
  
