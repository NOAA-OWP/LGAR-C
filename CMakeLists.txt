cmake_minimum_required(VERSION 3.10)

set(CMAKE_C_COMPILER $ENV{CC})
set(CMAKE_CXX_COMPILER $ENV{CXX})

string(ASCII 27 Esc)
set(ColourReset "${Esc}[m")
set(Red         "${Esc}[31m")
# module setup options

#option(BASE "BASE" ON)


# set the project name
project(lgarc)

set(CMAKE_BUILD_TYPE Debug)
IF(CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug build.")
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)
message(CMAKE_CXX_COMPILER " ${CMAKE_CXX_COMPILER}")
message(CMAKE_C_COMPILER " ${CMAKE_C_COMPILER}")
message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")

# -I option adds path to find header files
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I/usr/X11/include")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/X11/include")

# add the executable
message("CXXLIBS" ${CMAKE_CXX_FLAGS})

#add_executable(xlgar ./src/main.c ./src/lgar.c ./src/soil_funcs.c ./src/linked_list.c ./src/mem_funcs.c ./src/util_funcs.c)

#add_executable(xlgar ./src/main.cxx ./src/lgar.cxx ./src/soil_funcs.cxx ./src/linked_list.cxx ./src/mem_funcs.cxx ./src/util_funcs.cxx)

add_executable(xlgar ./src/bmi_main_lgar.cxx ./src/bmi_lgar.cxx ./src/lgar.cxx ./src/soil_funcs.cxx ./src/linked_list.cxx ./src/mem_funcs.cxx ./src/util_funcs.cxx)

#add_library(cfelib ./smc_coupler/src/bmi_soil_moisture_profile.cxx ./smc_coupler/src/soil_moisture_profile.cxx ./smc_coupler/include/bmi_soil_moisture_profile.hxx ./smc_coupler/include/soil_moisture_profile.hxx)


#target_link_libraries(xlgar LINK_PUBLIC xlgarlib)
target_link_libraries(xlgar PRIVATE m)

