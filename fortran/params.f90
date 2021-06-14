module params

  implicit none

  !logical to favor reducing the memory overhead at the cost of additional allocations
  logical, parameter :: favor_memory=.false.

  !logical to track the total energy while debugging
  logical, parameter :: track_energy=.false.
  
  !number of threads when executing in parallel or using OpenBLAS
  integer, parameter :: max_num_threads=16
  
  !real precision parameter
  integer, parameter :: rp=selected_real_kind(10,64)
  !integer precision parameter
  integer, parameter :: ip=selected_int_kind(5)

  !unit number for the program log and unit offset for other units
  !use unit_log=6 to output to terminal, otherwise make unit_log greater than 10
  !offset_units should be greater than unit_log
  integer, parameter :: unit_log=6, offset_units=20

contains

end module params
