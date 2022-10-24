cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_C_COMPILER $ENV{CC})
set(CMAKE_CXX_COMPILER $ENV{CXX})

string(ASCII 27 Esc)
set(ColourReset "${Esc}[m")
set(Red         "${Esc}[31m")
# module setup options

option(NGEN "NGEN" OFF)
option(STANDALONE "STANDALONE" OFF)
option(UNITTEST "UNITTEST" OFF)

if(NGEN)
  message("ngen framework build!")
  add_definitions(-DNGEN)
endif()
message("${NGEN} ${STANDALONE}")

if(STANDALONE)
 set(exe_name "lasam_standalone")
 message("Standalone build!")
endif()

if(UNITTEST)
 set(exe_name "lasam_unitest")
 message("Unittest build!")
endif()

# set the project name
project(lasambmi VERSION 1.0.0 DESCRIPTION "OWP LASAM BMI Module Shared Library")
#project(lgarc)

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
#message("CXXLIBS" ${CMAKE_CXX_FLAGS})


IF((NOT ${NGEN}) AND (NOT ${STANDALONE}))
  message("Pseudo-framework build")
ENDIF((NOT ${NGEN}) AND (NOT ${STANDALONE}))

if(STANDALONE)
  add_executable(${exe_name} ./src/bmi_main_lgar.cxx ./src/bmi_lgar.cxx ./src/lgar.cxx ./src/soil_funcs.cxx ./src/linked_list.cxx ./src/mem_funcs.cxx ./src/util_funcs.cxx ./src/aet.cxx ./giuh/giuh.h ./giuh/giuh.c)
  target_link_libraries(${exe_name} PRIVATE m)
elseif(UNITTEST)
  add_executable(${exe_name} ./test/main_unit_test_bmi.cxx ./src/bmi_lgar.cxx ./src/lgar.cxx ./src/soil_funcs.cxx ./src/linked_list.cxx ./src/mem_funcs.cxx ./src/util_funcs.cxx ./src/aet.cxx ./giuh/giuh.h ./giuh/giuh.c)
  target_link_libraries(${exe_name} PRIVATE m)
endif()


# ngen LASAM (add shared library)
set(LASAM_LIB_NAME_CMAKE lasambmi)
set(LASAM_LIB_DESC_CMAKE "OWP LASAM BMI Module Shared Library")

# Make sure these are compiled with this directive
add_compile_definitions(BMI_ACTIVE)

#add_library(lasambmi SHARED src/bmi_lgar.cxx src/lgar.cxx ./src/soil_funcs.cxx ./src/linked_list.cxx ./src/mem_funcs.cxx ./src/util_funcs.cxx ./src/aet.cxx ./giuh/giuh.c include/all.hxx ./giuh/giuh.h)
#target_include_directories(lasambmi PRIVATE include)

if(WIN32)
  add_library(lasambmi SHARED src/bmi_lgar.cxx src/lgar.cxx ./src/soil_funcs.cxx ./src/linked_list.cxx ./src/mem_funcs.cxx ./src/util_funcs.cxx ./src/aet.cxx ./giuh/giuh.c include/all.hxx ./giuh/giuh.h)
else()
   add_library(lasambmi SHARED src/bmi_lgar.cxx src/lgar.cxx ./src/soil_funcs.cxx ./src/linked_list.cxx ./src/mem_funcs.cxx ./src/util_funcs.cxx ./src/aet.cxx ./giuh/giuh.c include/all.hxx ./giuh/giuh.h)
endif()

target_include_directories(lasambmi PRIVATE include)

set_target_properties(lasambmi PROPERTIES VERSION ${PROJECT_VERSION})

set_target_properties(lasambmi PROPERTIES PUBLIC_HEADER ./include/bmi_lgar.hxx)

include(GNUInstallDirs)

install(TARGETS lasambmi
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

configure_file(lasambmi.pc.in lasambmi.pc @ONLY)

install(FILES ${CMAKE_BINARY_DIR}/lasambmi.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)


unset(STANDALONE CACHE)
unset(UNITTEST CACHE)
unset(NGEN CACHE)
