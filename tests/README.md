## Tests
#### - Unit test: Checks basic BMI functionality and runs model for one timestep (1 hour) to compare results against a benchmark.
#### - Synthetic tests: There are two synthetic examples for testing and to demonstrate the usage of the model. The examples use different soil hydraulic properties and different precipitation intensities without PET to simulate precipitation partitioning into infiltration and runoff. These simulate precipitation partitionind due to a short rainfall pulse that infiltrates entirely, followed by a longer precipitation pulse that generates runoff.
#### - Tests using observed forcing and soil hydraulic data: There are two examples using forcing and soil hydraulic data from USDA SCAN sites, near Phillipsburg, KS and Bushland, TX. Note that the outputs for these USDA examples are included in the tests folder but these are run from the 'LGAR-C' folder. The results of these simulations will be different than in the LGAR paper using these sites, as the forcing data in this repo have been resampled differently, and in the Phillipsburg example in this repo, the maximum ponded head is 2 cm, rather than 0 cm (as in the paper).


#### Build:
  - git clone https://github.com/NOAA-OWP/LGAR-C
  - cd LGAR-C && mkdir build && cd build
  - cmake ../ -DUNITTEST:BOOL=ON (for unittest)
  - cmake ../ -DSTANDALONE:BOOL=ON (for synthetic test)
  - make && cd ..

#### Run:
  - cd test
  - run `./run_unittest.sh` (for unittest)
  - run `./run_synthetic.sh OPTION` (for synthetic test; OPTION = 1 or 3 - these numbers correspond to different examples)

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
  $\textcolor{green}{\text{| All BMI Tests passed: Yes} }$ \
  $\textcolor{green}{\text{| Infiltration (mm) : (benchmark vs computed) | 1.896 vs 1.896} }$ \
  $\textcolor{green}{\text{| PET (mm): (benchmark vs computed) | 0.104 vs 0.104} }$ \
  $\textcolor{green}{\text{| AET (mm): (benchmark vs computed) | 0.0109182 vs 0.0109182} }$ \
  $\textcolor{green}{\text{| ************************************************************} }$


  Finally, here are the simulation results for the two synthetic cases and the two USDA SCAN cases:

  synth_1 should yield:

*********************************************************
-------------------- Simulation Summary -----------------
------------------------ Mass balance -------------------
Initial water in soil    =  24.0425193926 cm
Total precipitation      =  12.5000000000 cm
Total infiltration       =   4.4572080075 cm
Final water in soil      =  28.4997274000 cm
Surface ponded water     =   0.0000000000 cm
Surface runoff           =   8.0427919925 cm
GIUH runoff              =   0.0000000000 cm
Total percolation        =   0.0000000000 cm
Total AET                =   0.0000000000 cm
Total PET                =   0.0000000000 cm
Total discharge (Q)      =   0.0000000000 cm
Global balance           =   1.183764e-11 cm
Time                     =   0.01672 sec


synth_3 should yield:

*********************************************************
-------------------- Simulation Summary -----------------
------------------------ Mass balance -------------------
Initial water in soil    =  17.8650130977 cm
Total precipitation      =  31.2500000000 cm
Total infiltration       =  11.9438494670 cm
Final water in soil      =  29.8088625647 cm
Surface ponded water     =   0.0000000000 cm
Surface runoff           =  19.3061505330 cm
GIUH runoff              =   0.0000000000 cm
Total percolation        =   0.0000000000 cm
Total AET                =   0.0000000000 cm
Total PET                =   0.0000000000 cm
Total discharge (Q)      =   0.0000000000 cm
Global balance           =   3.394263e-11 cm
Time                     =   0.01874 sec


Phillipsburg should yield:

*********************************************************
-------------------- Simulation Summary -----------------
------------------------ Mass balance -------------------
Initial water in soil    =  45.1158503564 cm
Total precipitation      =  99.0854000000 cm
Total infiltration       =  88.6235782151 cm
Final water in soil      =  45.8124751752 cm
Surface ponded water     =   0.0000000000 cm
Surface runoff           =  10.4618217849 cm
GIUH runoff              =  10.4618217849 cm
Total percolation        =   0.0000000000 cm
Total AET                =  87.9269533606 cm
Total PET                = 147.9532940266 cm
Total discharge (Q)      =  10.4618217849 cm
Global balance           =   3.567497e-08 cm
Time                     =   5.374 sec


Bushland should yield:

*********************************************************
-------------------- Simulation Summary -----------------
------------------------ Mass balance -------------------
Initial water in soil    =  47.9782628218 cm
Total precipitation      =  25.9842000000 cm
Total infiltration       =  19.3791839646 cm
Final water in soil      =  42.8545804034 cm
Surface ponded water     =   0.0000000000 cm
Surface runoff           =   6.6050160354 cm
GIUH runoff              =   6.6050160354 cm
Total percolation        =   0.0000000000 cm
Total AET                =  24.5028663600 cm
Total PET                = 183.5308312160 cm
Total discharge (Q)      =   6.6050160354 cm
Global balance           =   2.305508e-08 cm
Time                     =   2.161 sec
