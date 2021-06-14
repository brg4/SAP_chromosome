# SAP_chromosome

Author: Benjamin R. Gilbert

Email: brg4@illinois.edu

### Directories ###

 fortran/ - contains source code for primary fortran program
 
 python/ - contains Jupyter notebooks used to pre-process inputs and post-process outputs of the fortran program
 
 example/ - contains an example set of input files

### Installation ###

 This program has been tested and built using gfortran-v9.3.0, packaged with the GCC-v9.3.0.

 You will need to specify the location of your installed fortran compiler using the variable "F90" in the Makefile.

 The program may optionally use an external BLAS library. Enable this by compiling the alternative program listed in the Makefile. The alternative program has been tested and built using OpenBLAS-v0.3.9.

 1) change directories to the fortran directory
 2) make clean
 3) make ring_generate.exe (or ring_generate_OBLAS.exe)

 Prior to building, the file params.f90 can be modified to tune the memory usage, change the maximum number of threads, specify the precision of real and integer variables, and the log file.

 The program can be run serially or in parallel using OpenMP. These options are specified in the input file.

### Execution ###

 The program is run by calling the executable "ring_generate.exe" on the command line. There are three optional command line arguments that can be included, all three are shown in the example below, and they can be included in any order.

 "./ring_generate.exe input=/some_input_directory/input_file.inp log=/some_log_directory/log_file.log num_threads=X"

  1) "input=(stuff here)" - specify the input file for the program
  2) "log=(stuff here)" - specify the log file for the program
  3) "num_threads=(a number)" - specify the number of threads to be used

  To run the example with two threads, change directories to the fortran directory, and enter the following on the command line

     "./ring_generate.exe input=../example/input/input_file.inp num_threads=2"

### Input Files ###

 There are four input files:

  1) main input file - ex. "input_file.inp" - all other input files are specified within this
   2) final species input file - ex. "final_species.dat"
   3) potential input file for growth - ex. "pot_grow.dat"
   4) potential input file for movement - ex. "pot_move.dat"

 Examples of all of these are provided within the example directory. Also within the example directory there is a file detailing the format of the input files.

### Output Files ###

 An output directory and output label are specified within the input file. All output files are written within the output directory and they all include the output label as a prefix.

 There are five output files:

  1) "(output_label).dat" - monomer indices and coordinates on integer lattice after successful termination
  2) "(output_label)_obst.dat" - placed obstacle indices, types, and coordinates on integer lattice after successful termination
  3) "(output_label)_species.dat" - monomer indices, species, and coordinates on integer lattice after successful termination
  4) "(output_label)_grow.dat" - monomer indices and coordinates on integer lattice at end of growth
  5) "(output_label)_species_grow.dat" - monomer indices and coordinates on integer lattice at end of growth

 In addition to the output files, there are intermediate restart files written after each error-checking routine is completed and a terminal restart file written after successful termination. The restart file will indicate if it is an intermediate restart file or a terminal restart file. All restart files match the format used by the input files, but will require minor modifications before they can be used, this is to prevent external scripts from creating infinite loops.

 Restart files are named as "(output_label)_restart.inp".

 If more than one replicate is being run, the output label is altered to indicate the individual replicates.

  output_label <= output_label+"_repXXXXX"

 For example, the output label for replicate 5 of a system with output_label="test" would become the following

  output_label="test_rep00005"

 and all output files and restart files for replicate 5 would use this output label.