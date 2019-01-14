cmake_minimum_required(VERSION 3.10.0)
project(ANA VERSION 1.1.0 LANGUAGES CXX)
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")

###################################################################
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()
###################################################################

# set(ANA_CGAL_DIR "" CACHE STRING "CGAL installation folder.")

add_executable(ANA "")
target_sources(ANA
    PRIVATE
        src/main.cpp
        src/ConvexHull.cpp
        src/main.cpp
        src/MD.cpp
        src/Modes.cpp
        src/NDD.cpp
        src/NDDUtils.cpp
        src/Primitives.cpp
        src/ProgramOptions.cpp
        src/Read.cpp
        src/Static.cpp
        src/Utils.cpp
        src/Write.cpp

        include/ANA/ConvexHull.hpp
        include/ANA/Includes.hpp
        include/ANA/MD.hpp
        include/ANA/Modes.hpp
        include/ANA/NDD.hpp
        include/ANA/NDDUtils.hpp
        include/ANA/Options.hpp
        include/ANA/pch.sh
        include/ANA/Primitives.hpp
        include/ANA/ProgramOptions.hpp
        include/ANA/Read.hpp
        include/ANA/Static.hpp
        include/ANA/Utils.hpp
        include/ANA/Write.hpp
)

message("${Magenta}Trying to automatically locate CGAL, gmp and mpfr libraries.${ColourReset}")
find_package(CGAL CONFIG REQUIRED)
target_link_libraries(ANA PRIVATE CGAL)
target_link_libraries(ANA PRIVATE gmp)
target_link_libraries(ANA PRIVATE mpfr)

message("${Magenta}Trying to automatically locate Boost.${ColourReset}")
find_package(Boost REQUIRED COMPONENTS program_options)
target_link_libraries(ANA PRIVATE boost_program_options)

message("${Magenta}Building own chemfiles library.${ColourReset}")
add_subdirectory(chemfiles EXCLUDE_FROM_ALL)
target_link_libraries(ANA PRIVATE chemfiles)

message("${Magenta}Using chemfiles's fmt header-only library.${ColourReset}")
include_directories(build/external/chemfiles/external/fmt)

target_include_directories(ANA
    PRIVATE
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

message("${Magenta}Including cotire for pre-compiled header.${ColourReset}")
include(cotire)
cotire(ANA)

target_compile_options(ANA PRIVATE -Wall -Wextra -Wpedantic)
target_compile_features(ANA PRIVATE cxx_std_17)
