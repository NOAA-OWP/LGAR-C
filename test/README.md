# Unit test for Lumped Arid/Semi-arid Model (LASAM)
Build instructions:
- mkdir build && cd build
- cmake ../ -DUNITTEST:BOOL=ON
- make

Running instructions:
- cd test
- run `./run_unittest.sh`

Several checks are performed:
1. Check names, number, and memory allocation of BMI input/output variables
2. Check units of BMI input/output variables 
3. Check number of layers and number of wetting fronts against the benchmark
4. Check BMI Grid function (e.g., grid_id, grid_size, etc.)
5. Test `GetVar*` methods for the BMI input variables and compare against initial (or prescribed) data; if failed, will throw an error
6. Loop over the input variables, use `Set*` and `Get*` methods to verify `Get*` return the same data set by `Set*`
7. Using `Update` method, advance the model to get updated depths and soil moisture of the wetting fronts. Compare against the benchmark values.
