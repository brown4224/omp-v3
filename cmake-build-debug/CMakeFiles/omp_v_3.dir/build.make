# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.8

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/sean/Desktop/clion/bin/cmake/bin/cmake

# The command to remove a file.
RM = /home/sean/Desktop/clion/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/sean/Desktop/school/parralel-systems/omp-v-3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/sean/Desktop/school/parralel-systems/omp-v-3/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/omp_v_3.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/omp_v_3.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/omp_v_3.dir/flags.make

CMakeFiles/omp_v_3.dir/main.cpp.o: CMakeFiles/omp_v_3.dir/flags.make
CMakeFiles/omp_v_3.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sean/Desktop/school/parralel-systems/omp-v-3/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/omp_v_3.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/omp_v_3.dir/main.cpp.o -c /home/sean/Desktop/school/parralel-systems/omp-v-3/main.cpp

CMakeFiles/omp_v_3.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/omp_v_3.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sean/Desktop/school/parralel-systems/omp-v-3/main.cpp > CMakeFiles/omp_v_3.dir/main.cpp.i

CMakeFiles/omp_v_3.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/omp_v_3.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sean/Desktop/school/parralel-systems/omp-v-3/main.cpp -o CMakeFiles/omp_v_3.dir/main.cpp.s

CMakeFiles/omp_v_3.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/omp_v_3.dir/main.cpp.o.requires

CMakeFiles/omp_v_3.dir/main.cpp.o.provides: CMakeFiles/omp_v_3.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/omp_v_3.dir/build.make CMakeFiles/omp_v_3.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/omp_v_3.dir/main.cpp.o.provides

CMakeFiles/omp_v_3.dir/main.cpp.o.provides.build: CMakeFiles/omp_v_3.dir/main.cpp.o


# Object files for target omp_v_3
omp_v_3_OBJECTS = \
"CMakeFiles/omp_v_3.dir/main.cpp.o"

# External object files for target omp_v_3
omp_v_3_EXTERNAL_OBJECTS =

omp_v_3: CMakeFiles/omp_v_3.dir/main.cpp.o
omp_v_3: CMakeFiles/omp_v_3.dir/build.make
omp_v_3: CMakeFiles/omp_v_3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/sean/Desktop/school/parralel-systems/omp-v-3/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable omp_v_3"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/omp_v_3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/omp_v_3.dir/build: omp_v_3

.PHONY : CMakeFiles/omp_v_3.dir/build

CMakeFiles/omp_v_3.dir/requires: CMakeFiles/omp_v_3.dir/main.cpp.o.requires

.PHONY : CMakeFiles/omp_v_3.dir/requires

CMakeFiles/omp_v_3.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/omp_v_3.dir/cmake_clean.cmake
.PHONY : CMakeFiles/omp_v_3.dir/clean

CMakeFiles/omp_v_3.dir/depend:
	cd /home/sean/Desktop/school/parralel-systems/omp-v-3/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/sean/Desktop/school/parralel-systems/omp-v-3 /home/sean/Desktop/school/parralel-systems/omp-v-3 /home/sean/Desktop/school/parralel-systems/omp-v-3/cmake-build-debug /home/sean/Desktop/school/parralel-systems/omp-v-3/cmake-build-debug /home/sean/Desktop/school/parralel-systems/omp-v-3/cmake-build-debug/CMakeFiles/omp_v_3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/omp_v_3.dir/depend

