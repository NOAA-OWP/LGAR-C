## Unit test for Lumped Arid/Semi-arid Model (LASAM)
#### Build:
  - mkdir build && cd build
  - cmake ../ -DUNITTEST:BOOL=ON
  - make

#### Run:
  - cd test
  - run `./run_unittest.sh`

#### Several checks are performed:
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
  $\textcolor{green}{\text{| Infiltration: (benchmark vs computed) | 0.173583 vs 0.173583} }$ \
  $\textcolor{green}{\text{| PET: (benchmark vs computed) | 0.0104167 vs 0.0104167} }$ \
  $\textcolor{green}{\text{| AET: (benchmark vs computed) | 0.000845446 vs 0.000845446} }$ \
  $\textcolor{green}{\text{| ************************************************************} }$
  
