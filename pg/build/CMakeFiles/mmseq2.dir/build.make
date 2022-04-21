# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/michal/bachelors-thesis/bioseqdb/pg

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/michal/bachelors-thesis/bioseqdb/pg/build

# Include any dependencies generated for this target.
include CMakeFiles/mmseq2.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/mmseq2.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/mmseq2.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mmseq2.dir/flags.make

CMakeFiles/mmseq2.dir/extension.cpp.o: CMakeFiles/mmseq2.dir/flags.make
CMakeFiles/mmseq2.dir/extension.cpp.o: ../extension.cpp
CMakeFiles/mmseq2.dir/extension.cpp.o: CMakeFiles/mmseq2.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/michal/bachelors-thesis/bioseqdb/pg/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/mmseq2.dir/extension.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/mmseq2.dir/extension.cpp.o -MF CMakeFiles/mmseq2.dir/extension.cpp.o.d -o CMakeFiles/mmseq2.dir/extension.cpp.o -c /home/michal/bachelors-thesis/bioseqdb/pg/extension.cpp

CMakeFiles/mmseq2.dir/extension.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mmseq2.dir/extension.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/michal/bachelors-thesis/bioseqdb/pg/extension.cpp > CMakeFiles/mmseq2.dir/extension.cpp.i

CMakeFiles/mmseq2.dir/extension.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mmseq2.dir/extension.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/michal/bachelors-thesis/bioseqdb/pg/extension.cpp -o CMakeFiles/mmseq2.dir/extension.cpp.s

# Object files for target mmseq2
mmseq2_OBJECTS = \
"CMakeFiles/mmseq2.dir/extension.cpp.o"

# External object files for target mmseq2
mmseq2_EXTERNAL_OBJECTS =

libmmseq2.so: CMakeFiles/mmseq2.dir/extension.cpp.o
libmmseq2.so: CMakeFiles/mmseq2.dir/build.make
libmmseq2.so: /usr/lib/x86_64-linux-gnu/libpq.so
libmmseq2.so: CMakeFiles/mmseq2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/michal/bachelors-thesis/bioseqdb/pg/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library libmmseq2.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mmseq2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mmseq2.dir/build: libmmseq2.so
.PHONY : CMakeFiles/mmseq2.dir/build

CMakeFiles/mmseq2.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mmseq2.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mmseq2.dir/clean

CMakeFiles/mmseq2.dir/depend:
	cd /home/michal/bachelors-thesis/bioseqdb/pg/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/michal/bachelors-thesis/bioseqdb/pg /home/michal/bachelors-thesis/bioseqdb/pg /home/michal/bachelors-thesis/bioseqdb/pg/build /home/michal/bachelors-thesis/bioseqdb/pg/build /home/michal/bachelors-thesis/bioseqdb/pg/build/CMakeFiles/mmseq2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mmseq2.dir/depend

