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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/joshwolper/Desktop/CIS563_SmokeBaseCode

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug

# Include any dependencies generated for this target.
include src/CMakeFiles/SMOKE.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/SMOKE.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/SMOKE.dir/flags.make

src/CMakeFiles/SMOKE.dir/main.cpp.o: src/CMakeFiles/SMOKE.dir/flags.make
src/CMakeFiles/SMOKE.dir/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/SMOKE.dir/main.cpp.o"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SMOKE.dir/main.cpp.o -c /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/main.cpp

src/CMakeFiles/SMOKE.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SMOKE.dir/main.cpp.i"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/main.cpp > CMakeFiles/SMOKE.dir/main.cpp.i

src/CMakeFiles/SMOKE.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SMOKE.dir/main.cpp.s"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/main.cpp -o CMakeFiles/SMOKE.dir/main.cpp.s

src/CMakeFiles/SMOKE.dir/main.cpp.o.requires:

.PHONY : src/CMakeFiles/SMOKE.dir/main.cpp.o.requires

src/CMakeFiles/SMOKE.dir/main.cpp.o.provides: src/CMakeFiles/SMOKE.dir/main.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/SMOKE.dir/build.make src/CMakeFiles/SMOKE.dir/main.cpp.o.provides.build
.PHONY : src/CMakeFiles/SMOKE.dir/main.cpp.o.provides

src/CMakeFiles/SMOKE.dir/main.cpp.o.provides.build: src/CMakeFiles/SMOKE.dir/main.cpp.o


src/CMakeFiles/SMOKE.dir/smoke_sim.cpp.o: src/CMakeFiles/SMOKE.dir/flags.make
src/CMakeFiles/SMOKE.dir/smoke_sim.cpp.o: ../src/smoke_sim.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/CMakeFiles/SMOKE.dir/smoke_sim.cpp.o"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SMOKE.dir/smoke_sim.cpp.o -c /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/smoke_sim.cpp

src/CMakeFiles/SMOKE.dir/smoke_sim.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SMOKE.dir/smoke_sim.cpp.i"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/smoke_sim.cpp > CMakeFiles/SMOKE.dir/smoke_sim.cpp.i

src/CMakeFiles/SMOKE.dir/smoke_sim.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SMOKE.dir/smoke_sim.cpp.s"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/smoke_sim.cpp -o CMakeFiles/SMOKE.dir/smoke_sim.cpp.s

src/CMakeFiles/SMOKE.dir/smoke_sim.cpp.o.requires:

.PHONY : src/CMakeFiles/SMOKE.dir/smoke_sim.cpp.o.requires

src/CMakeFiles/SMOKE.dir/smoke_sim.cpp.o.provides: src/CMakeFiles/SMOKE.dir/smoke_sim.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/SMOKE.dir/build.make src/CMakeFiles/SMOKE.dir/smoke_sim.cpp.o.provides.build
.PHONY : src/CMakeFiles/SMOKE.dir/smoke_sim.cpp.o.provides

src/CMakeFiles/SMOKE.dir/smoke_sim.cpp.o.provides.build: src/CMakeFiles/SMOKE.dir/smoke_sim.cpp.o


src/CMakeFiles/SMOKE.dir/mac_grid.cpp.o: src/CMakeFiles/SMOKE.dir/flags.make
src/CMakeFiles/SMOKE.dir/mac_grid.cpp.o: ../src/mac_grid.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/CMakeFiles/SMOKE.dir/mac_grid.cpp.o"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SMOKE.dir/mac_grid.cpp.o -c /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/mac_grid.cpp

src/CMakeFiles/SMOKE.dir/mac_grid.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SMOKE.dir/mac_grid.cpp.i"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/mac_grid.cpp > CMakeFiles/SMOKE.dir/mac_grid.cpp.i

src/CMakeFiles/SMOKE.dir/mac_grid.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SMOKE.dir/mac_grid.cpp.s"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/mac_grid.cpp -o CMakeFiles/SMOKE.dir/mac_grid.cpp.s

src/CMakeFiles/SMOKE.dir/mac_grid.cpp.o.requires:

.PHONY : src/CMakeFiles/SMOKE.dir/mac_grid.cpp.o.requires

src/CMakeFiles/SMOKE.dir/mac_grid.cpp.o.provides: src/CMakeFiles/SMOKE.dir/mac_grid.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/SMOKE.dir/build.make src/CMakeFiles/SMOKE.dir/mac_grid.cpp.o.provides.build
.PHONY : src/CMakeFiles/SMOKE.dir/mac_grid.cpp.o.provides

src/CMakeFiles/SMOKE.dir/mac_grid.cpp.o.provides.build: src/CMakeFiles/SMOKE.dir/mac_grid.cpp.o


src/CMakeFiles/SMOKE.dir/vec.cpp.o: src/CMakeFiles/SMOKE.dir/flags.make
src/CMakeFiles/SMOKE.dir/vec.cpp.o: ../src/vec.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/CMakeFiles/SMOKE.dir/vec.cpp.o"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SMOKE.dir/vec.cpp.o -c /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/vec.cpp

src/CMakeFiles/SMOKE.dir/vec.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SMOKE.dir/vec.cpp.i"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/vec.cpp > CMakeFiles/SMOKE.dir/vec.cpp.i

src/CMakeFiles/SMOKE.dir/vec.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SMOKE.dir/vec.cpp.s"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/vec.cpp -o CMakeFiles/SMOKE.dir/vec.cpp.s

src/CMakeFiles/SMOKE.dir/vec.cpp.o.requires:

.PHONY : src/CMakeFiles/SMOKE.dir/vec.cpp.o.requires

src/CMakeFiles/SMOKE.dir/vec.cpp.o.provides: src/CMakeFiles/SMOKE.dir/vec.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/SMOKE.dir/build.make src/CMakeFiles/SMOKE.dir/vec.cpp.o.provides.build
.PHONY : src/CMakeFiles/SMOKE.dir/vec.cpp.o.provides

src/CMakeFiles/SMOKE.dir/vec.cpp.o.provides.build: src/CMakeFiles/SMOKE.dir/vec.cpp.o


src/CMakeFiles/SMOKE.dir/grid_data.cpp.o: src/CMakeFiles/SMOKE.dir/flags.make
src/CMakeFiles/SMOKE.dir/grid_data.cpp.o: ../src/grid_data.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object src/CMakeFiles/SMOKE.dir/grid_data.cpp.o"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SMOKE.dir/grid_data.cpp.o -c /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/grid_data.cpp

src/CMakeFiles/SMOKE.dir/grid_data.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SMOKE.dir/grid_data.cpp.i"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/grid_data.cpp > CMakeFiles/SMOKE.dir/grid_data.cpp.i

src/CMakeFiles/SMOKE.dir/grid_data.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SMOKE.dir/grid_data.cpp.s"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/grid_data.cpp -o CMakeFiles/SMOKE.dir/grid_data.cpp.s

src/CMakeFiles/SMOKE.dir/grid_data.cpp.o.requires:

.PHONY : src/CMakeFiles/SMOKE.dir/grid_data.cpp.o.requires

src/CMakeFiles/SMOKE.dir/grid_data.cpp.o.provides: src/CMakeFiles/SMOKE.dir/grid_data.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/SMOKE.dir/build.make src/CMakeFiles/SMOKE.dir/grid_data.cpp.o.provides.build
.PHONY : src/CMakeFiles/SMOKE.dir/grid_data.cpp.o.provides

src/CMakeFiles/SMOKE.dir/grid_data.cpp.o.provides.build: src/CMakeFiles/SMOKE.dir/grid_data.cpp.o


src/CMakeFiles/SMOKE.dir/camera.cpp.o: src/CMakeFiles/SMOKE.dir/flags.make
src/CMakeFiles/SMOKE.dir/camera.cpp.o: ../src/camera.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object src/CMakeFiles/SMOKE.dir/camera.cpp.o"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SMOKE.dir/camera.cpp.o -c /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/camera.cpp

src/CMakeFiles/SMOKE.dir/camera.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SMOKE.dir/camera.cpp.i"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/camera.cpp > CMakeFiles/SMOKE.dir/camera.cpp.i

src/CMakeFiles/SMOKE.dir/camera.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SMOKE.dir/camera.cpp.s"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/camera.cpp -o CMakeFiles/SMOKE.dir/camera.cpp.s

src/CMakeFiles/SMOKE.dir/camera.cpp.o.requires:

.PHONY : src/CMakeFiles/SMOKE.dir/camera.cpp.o.requires

src/CMakeFiles/SMOKE.dir/camera.cpp.o.provides: src/CMakeFiles/SMOKE.dir/camera.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/SMOKE.dir/build.make src/CMakeFiles/SMOKE.dir/camera.cpp.o.provides.build
.PHONY : src/CMakeFiles/SMOKE.dir/camera.cpp.o.provides

src/CMakeFiles/SMOKE.dir/camera.cpp.o.provides.build: src/CMakeFiles/SMOKE.dir/camera.cpp.o


src/CMakeFiles/SMOKE.dir/fps.cpp.o: src/CMakeFiles/SMOKE.dir/flags.make
src/CMakeFiles/SMOKE.dir/fps.cpp.o: ../src/fps.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object src/CMakeFiles/SMOKE.dir/fps.cpp.o"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SMOKE.dir/fps.cpp.o -c /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/fps.cpp

src/CMakeFiles/SMOKE.dir/fps.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SMOKE.dir/fps.cpp.i"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/fps.cpp > CMakeFiles/SMOKE.dir/fps.cpp.i

src/CMakeFiles/SMOKE.dir/fps.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SMOKE.dir/fps.cpp.s"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/fps.cpp -o CMakeFiles/SMOKE.dir/fps.cpp.s

src/CMakeFiles/SMOKE.dir/fps.cpp.o.requires:

.PHONY : src/CMakeFiles/SMOKE.dir/fps.cpp.o.requires

src/CMakeFiles/SMOKE.dir/fps.cpp.o.provides: src/CMakeFiles/SMOKE.dir/fps.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/SMOKE.dir/build.make src/CMakeFiles/SMOKE.dir/fps.cpp.o.provides.build
.PHONY : src/CMakeFiles/SMOKE.dir/fps.cpp.o.provides

src/CMakeFiles/SMOKE.dir/fps.cpp.o.provides.build: src/CMakeFiles/SMOKE.dir/fps.cpp.o


src/CMakeFiles/SMOKE.dir/constants.cpp.o: src/CMakeFiles/SMOKE.dir/flags.make
src/CMakeFiles/SMOKE.dir/constants.cpp.o: ../src/constants.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object src/CMakeFiles/SMOKE.dir/constants.cpp.o"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SMOKE.dir/constants.cpp.o -c /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/constants.cpp

src/CMakeFiles/SMOKE.dir/constants.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SMOKE.dir/constants.cpp.i"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/constants.cpp > CMakeFiles/SMOKE.dir/constants.cpp.i

src/CMakeFiles/SMOKE.dir/constants.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SMOKE.dir/constants.cpp.s"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/constants.cpp -o CMakeFiles/SMOKE.dir/constants.cpp.s

src/CMakeFiles/SMOKE.dir/constants.cpp.o.requires:

.PHONY : src/CMakeFiles/SMOKE.dir/constants.cpp.o.requires

src/CMakeFiles/SMOKE.dir/constants.cpp.o.provides: src/CMakeFiles/SMOKE.dir/constants.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/SMOKE.dir/build.make src/CMakeFiles/SMOKE.dir/constants.cpp.o.provides.build
.PHONY : src/CMakeFiles/SMOKE.dir/constants.cpp.o.provides

src/CMakeFiles/SMOKE.dir/constants.cpp.o.provides.build: src/CMakeFiles/SMOKE.dir/constants.cpp.o


src/CMakeFiles/SMOKE.dir/basic_math.cpp.o: src/CMakeFiles/SMOKE.dir/flags.make
src/CMakeFiles/SMOKE.dir/basic_math.cpp.o: ../src/basic_math.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object src/CMakeFiles/SMOKE.dir/basic_math.cpp.o"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SMOKE.dir/basic_math.cpp.o -c /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/basic_math.cpp

src/CMakeFiles/SMOKE.dir/basic_math.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SMOKE.dir/basic_math.cpp.i"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/basic_math.cpp > CMakeFiles/SMOKE.dir/basic_math.cpp.i

src/CMakeFiles/SMOKE.dir/basic_math.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SMOKE.dir/basic_math.cpp.s"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src/basic_math.cpp -o CMakeFiles/SMOKE.dir/basic_math.cpp.s

src/CMakeFiles/SMOKE.dir/basic_math.cpp.o.requires:

.PHONY : src/CMakeFiles/SMOKE.dir/basic_math.cpp.o.requires

src/CMakeFiles/SMOKE.dir/basic_math.cpp.o.provides: src/CMakeFiles/SMOKE.dir/basic_math.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/SMOKE.dir/build.make src/CMakeFiles/SMOKE.dir/basic_math.cpp.o.provides.build
.PHONY : src/CMakeFiles/SMOKE.dir/basic_math.cpp.o.provides

src/CMakeFiles/SMOKE.dir/basic_math.cpp.o.provides.build: src/CMakeFiles/SMOKE.dir/basic_math.cpp.o


# Object files for target SMOKE
SMOKE_OBJECTS = \
"CMakeFiles/SMOKE.dir/main.cpp.o" \
"CMakeFiles/SMOKE.dir/smoke_sim.cpp.o" \
"CMakeFiles/SMOKE.dir/mac_grid.cpp.o" \
"CMakeFiles/SMOKE.dir/vec.cpp.o" \
"CMakeFiles/SMOKE.dir/grid_data.cpp.o" \
"CMakeFiles/SMOKE.dir/camera.cpp.o" \
"CMakeFiles/SMOKE.dir/fps.cpp.o" \
"CMakeFiles/SMOKE.dir/constants.cpp.o" \
"CMakeFiles/SMOKE.dir/basic_math.cpp.o"

# External object files for target SMOKE
SMOKE_EXTERNAL_OBJECTS =

../src/SMOKE_debug: src/CMakeFiles/SMOKE.dir/main.cpp.o
../src/SMOKE_debug: src/CMakeFiles/SMOKE.dir/smoke_sim.cpp.o
../src/SMOKE_debug: src/CMakeFiles/SMOKE.dir/mac_grid.cpp.o
../src/SMOKE_debug: src/CMakeFiles/SMOKE.dir/vec.cpp.o
../src/SMOKE_debug: src/CMakeFiles/SMOKE.dir/grid_data.cpp.o
../src/SMOKE_debug: src/CMakeFiles/SMOKE.dir/camera.cpp.o
../src/SMOKE_debug: src/CMakeFiles/SMOKE.dir/fps.cpp.o
../src/SMOKE_debug: src/CMakeFiles/SMOKE.dir/constants.cpp.o
../src/SMOKE_debug: src/CMakeFiles/SMOKE.dir/basic_math.cpp.o
../src/SMOKE_debug: src/CMakeFiles/SMOKE.dir/build.make
../src/SMOKE_debug: partio-build/lib/libpartio.a
../src/SMOKE_debug: src/CMakeFiles/SMOKE.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Linking CXX executable ../../src/SMOKE_debug"
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SMOKE.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/SMOKE.dir/build: ../src/SMOKE_debug

.PHONY : src/CMakeFiles/SMOKE.dir/build

src/CMakeFiles/SMOKE.dir/requires: src/CMakeFiles/SMOKE.dir/main.cpp.o.requires
src/CMakeFiles/SMOKE.dir/requires: src/CMakeFiles/SMOKE.dir/smoke_sim.cpp.o.requires
src/CMakeFiles/SMOKE.dir/requires: src/CMakeFiles/SMOKE.dir/mac_grid.cpp.o.requires
src/CMakeFiles/SMOKE.dir/requires: src/CMakeFiles/SMOKE.dir/vec.cpp.o.requires
src/CMakeFiles/SMOKE.dir/requires: src/CMakeFiles/SMOKE.dir/grid_data.cpp.o.requires
src/CMakeFiles/SMOKE.dir/requires: src/CMakeFiles/SMOKE.dir/camera.cpp.o.requires
src/CMakeFiles/SMOKE.dir/requires: src/CMakeFiles/SMOKE.dir/fps.cpp.o.requires
src/CMakeFiles/SMOKE.dir/requires: src/CMakeFiles/SMOKE.dir/constants.cpp.o.requires
src/CMakeFiles/SMOKE.dir/requires: src/CMakeFiles/SMOKE.dir/basic_math.cpp.o.requires

.PHONY : src/CMakeFiles/SMOKE.dir/requires

src/CMakeFiles/SMOKE.dir/clean:
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src && $(CMAKE_COMMAND) -P CMakeFiles/SMOKE.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/SMOKE.dir/clean

src/CMakeFiles/SMOKE.dir/depend:
	cd /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/joshwolper/Desktop/CIS563_SmokeBaseCode /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/src /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src /Users/joshwolper/Desktop/CIS563_SmokeBaseCode/cmake-build-debug/src/CMakeFiles/SMOKE.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/SMOKE.dir/depend
