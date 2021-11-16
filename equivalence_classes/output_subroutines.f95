
!
!   output_subroutines.f95
!
!   Generates all minimum size instances of the pallet loading problem for
!   less than N+1 boxes (by area).
!
!   Output subroutines.
!
!   For details on the algorithm, see 
!   "An algorithm for identifying equivalence classes of the pallet loading problem"
!
!   Ed Herbert, 12 Nov 2021
!
!


  !
  !   subroutine output_msi
  !

  subroutine output_msi(ix, iy, ia, ib)

    implicit none
    integer(8), intent(in) :: ix, iy, ia, ib
    integer(8) :: i

    Nobs = Nobs + 1

    if (output_problems .eq. 1) write (8, fmt = out_format) ix, iy, ia, ib

    if (output_summary .eq. 1) then
      ! integer division
      i = (ix * iy) / (ia * ib)
      summary(i, 1) = summary(i, 1) + 1
      if (ix .ge. summary(i, 2)) then
	summary(i, 2) = max(summary(i, 2), ix)
	maxx(i, 1) = ix
	maxx(i, 2) = iy
	maxx(i, 3) = ia
	maxx(i, 4) = ib
      end if
      if (iy .ge. summary(i, 3)) then
	summary(i, 3) = max(summary(i, 3), iy)
	maxy(i, 1) = ix
	maxy(i, 2) = iy
	maxy(i, 3) = ia
	maxy(i, 4) = ib
      end if
      if (ia .ge. summary(i, 4)) then
	summary(i, 4) = max(summary(i, 4), ia)
	maxa(i, 1) = ix
	maxa(i, 2) = iy
	maxa(i, 3) = ia
	maxa(i, 4) = ib
      end if
      if (ib .ge. summary(i, 5)) then
	summary(i, 5) = max(summary(i, 5), ib)
	maxb(i, 1) = ix
	maxb(i, 2) = iy
	maxb(i, 3) = ia
	maxb(i, 4) = ib
      end if

      if (ix .lt. summary(i, 6)) summary(i, 6) = ix
      if (iy .lt. summary(i, 7)) summary(i, 7) = iy
      if (ia .lt. summary(i, 8)) summary(i, 8) = ia
      if (ib .lt. summary(i, 9)) summary(i, 9) = ib

    end if

  end subroutine output_msi




  subroutine output_preamble(N)
  
    implicit none
    integer(8) :: N, itime
    character(len=30) :: date
   
    if (output_problems .eq. 1) open (unit = 8, file = "output.dat", status = "replace", action="write")
    if (output_summary .eq. 1) then
      open (unit = 10, file = "summary.dat", status = "replace", action="write")
      summary(1:N, 1:5) = 0
      summary(1:N, 6:9) = huge(N)
      maxa = 0
      maxb = 0
      maxx = 0
      maxy = 0
    end if

    itime = time8()
    call ctime(itime, date)

    write (unit = 9, fmt = "(a10i5a4a30)") 'Starting  ', N, " at ", date

  end subroutine output_preamble




  subroutine output_final(N, start, finish)
  
    implicit none
    integer(8) :: N, itime, i, total 
    real :: start, finish
    character(len=30) :: end_date

    itime = time8()
    call ctime(itime, end_date)
    write (unit = 9, fmt = "(a10i5a4a30a8f12.3a9i20a8)") 'Completed ', N, " at ", end_date, &
      " Time = ", finish - start, " seconds ", Nobs, " classes"

    if (output_problems .eq. 1) close (unit = 8)
    if (output_summary .eq. 1) then
      total = 0
      write (10, *) " N      Cum Observations      Observations     Max(X)    Max(Y) Max(a) Max(b)"
      do i = 1, N
        total = total + summary(i, 1)
        write (10, fmt = "(i5i20i20i12i7i7i20i12i7i7i7i12i7i7i7i12i7i7i7i12i7i7i7i12i7i7i7)") i, total, summary(i, 1:9), &
	    maxx(i, 1:4), maxy(i, 1:4), maxa(i, 1:4), maxb(i, 1:4)
      end do
      close (unit = 10)

    end if

  end subroutine output_final



