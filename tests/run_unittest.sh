#!/bin/bash
#${CXX} -lm -Wall -O -g ../src/bmi_main_lgar.cxx ../src/bmi_lgar.cxx ../src/lgar.cxx ../src/soil_funcs.cxx ../src/linked_list.cxx ../src/mem_funcs.cxx ../src/util_funcs.cxx ../src/aet.cxx ../giuh/giuh.h ../giuh/giuh.c -o run_lasam
../build/lasam_unitest configs/unittest.txt
