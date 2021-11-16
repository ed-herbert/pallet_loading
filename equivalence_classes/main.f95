
!
!   main.f95
!
!   Generates all minimum size instances of the pallet loading problem for
!   less than N+1 boxes (by area).
!
!   For details on the algorithm, see 
!   "An algorithm for identifying equivalence classes of the pallet loading problem"
!
!   Ed Herbert, 12 Nov 2021
!
!


Module generate_msi

  integer, parameter :: N_MAX = 600
  integer, parameter :: X_MAX = N_MAX*(N_MAX+1)
  integer, parameter :: X_MAX2 = 1 + X_MAX / 64 

  character(len = 15) :: out_format = "(i7xi4xi4xi4)"

  integer(8) :: boxes(2,2)
  integer(8) :: ixmin, ixmax, iymin, iymax, iamax, iamin, box_area
  integer(8) :: ix2, ixmin2, iymax2
  integer(8) :: ial, iau, ibl, ibu

  integer(8), dimension(X_MAX) :: izl
  integer(8) :: N_izl, min_area, max_area, Nobs

  integer(8), dimension(X_MAX2) :: is_tight, is_slack, is_tight_low, is_tight_high

  ! output variables
  integer(8) :: output_problems, output_summary
  integer(8) :: summary(N_MAX, 9)
  integer(8) :: maxa(N_MAX, 5)
  integer(8) :: maxb(N_MAX, 5)
  integer(8) :: maxx(N_MAX, 5)
  integer(8) :: maxy(N_MAX, 5)


  contains

  include 'generate_all_msi.f95'

  include 'output_subroutines.f95'

  include 'misc_functions.f95'

end module generate_msi




program msi_algorithm

  use generate_msi

  implicit none
  integer(8) :: N
  real :: start, finish
  character(len=30) :: start_date, end_date


  print *, "Generate all pallet loading problem instances"
  print *, "Instances file output.dat, summary file summary.dat, timings file timings.dat"
  print *, "Enter N (Max", N_MAX, "), Output all instances? (1/0), Output summary stats (1/0)"
  read *, N, output_problems, output_summary

  if (N .gt. N_MAX) then
    print *, "N greater than N_MAX. Either recompile or supply smaller N"
    print *, "N_MAX = ", N_MAX
    stop
  end if

  open (unit = 9, file = "timings.dat", status = "unknown", action="write", position="append")


  call output_preamble(N)

  call cpu_time(start)
  call generate_all_msi(N)
  call cpu_time(finish)

  call output_final(N, start, finish)


end program msi_algorithm


