# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/lrisser/Softwares/TV2016

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/lrisser/Softwares/TV2016/Build

# Include any dependencies generated for this target.
include lib/src/CMakeFiles/SciCalcPack.dir/depend.make

# Include the progress variables for this target.
include lib/src/CMakeFiles/SciCalcPack.dir/progress.make

# Include the compile flags for this target's objects.
include lib/src/CMakeFiles/SciCalcPack.dir/flags.make

lib/src/CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.o: lib/src/CMakeFiles/SciCalcPack.dir/flags.make
lib/src/CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.o: ../lib/src/SciCalcPack.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/lrisser/Softwares/TV2016/Build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object lib/src/CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.o"
	cd /home/lrisser/Softwares/TV2016/Build/lib/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.o -c /home/lrisser/Softwares/TV2016/lib/src/SciCalcPack.cc

lib/src/CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.i"
	cd /home/lrisser/Softwares/TV2016/Build/lib/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/lrisser/Softwares/TV2016/lib/src/SciCalcPack.cc > CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.i

lib/src/CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.s"
	cd /home/lrisser/Softwares/TV2016/Build/lib/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/lrisser/Softwares/TV2016/lib/src/SciCalcPack.cc -o CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.s

lib/src/CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.o.requires:
.PHONY : lib/src/CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.o.requires

lib/src/CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.o.provides: lib/src/CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.o.requires
	$(MAKE) -f lib/src/CMakeFiles/SciCalcPack.dir/build.make lib/src/CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.o.provides.build
.PHONY : lib/src/CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.o.provides

lib/src/CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.o.provides.build: lib/src/CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.o

# Object files for target SciCalcPack
SciCalcPack_OBJECTS = \
"CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.o"

# External object files for target SciCalcPack
SciCalcPack_EXTERNAL_OBJECTS =

lib/src/libSciCalcPack.a: lib/src/CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.o
lib/src/libSciCalcPack.a: lib/src/CMakeFiles/SciCalcPack.dir/build.make
lib/src/libSciCalcPack.a: lib/src/CMakeFiles/SciCalcPack.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libSciCalcPack.a"
	cd /home/lrisser/Softwares/TV2016/Build/lib/src && $(CMAKE_COMMAND) -P CMakeFiles/SciCalcPack.dir/cmake_clean_target.cmake
	cd /home/lrisser/Softwares/TV2016/Build/lib/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SciCalcPack.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/src/CMakeFiles/SciCalcPack.dir/build: lib/src/libSciCalcPack.a
.PHONY : lib/src/CMakeFiles/SciCalcPack.dir/build

lib/src/CMakeFiles/SciCalcPack.dir/requires: lib/src/CMakeFiles/SciCalcPack.dir/SciCalcPack.cc.o.requires
.PHONY : lib/src/CMakeFiles/SciCalcPack.dir/requires

lib/src/CMakeFiles/SciCalcPack.dir/clean:
	cd /home/lrisser/Softwares/TV2016/Build/lib/src && $(CMAKE_COMMAND) -P CMakeFiles/SciCalcPack.dir/cmake_clean.cmake
.PHONY : lib/src/CMakeFiles/SciCalcPack.dir/clean

lib/src/CMakeFiles/SciCalcPack.dir/depend:
	cd /home/lrisser/Softwares/TV2016/Build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lrisser/Softwares/TV2016 /home/lrisser/Softwares/TV2016/lib/src /home/lrisser/Softwares/TV2016/Build /home/lrisser/Softwares/TV2016/Build/lib/src /home/lrisser/Softwares/TV2016/Build/lib/src/CMakeFiles/SciCalcPack.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/src/CMakeFiles/SciCalcPack.dir/depend

