
!
!   generate_all_msi.f95
!
!   Generates all minimum size instances of the pallet loading problem for
!   less than N+1 boxes (by area).
!
!   Subroutines for doing the heavy lifting with Level III constraints.
!
!   For details on the algorithm, see 
!   "An algorithm for identifying equivalence classes of the pallet loading problem"
!
!   Ed Herbert, 12 Nov 2021
!
!





  !
  !  subroutine generate_all_msi
  !
  !  Generates all MSI classes for .le. N boxes by area bound
  !

  subroutine generate_all_msi(N)

    implicit none
    integer(8), intent(in) :: N
    integer(8) ia, ib, ix, iy
    integer(8) ix_index, iy_index, min_area, ipallet

    Nobs = 0

    do ib = 1, N

      !
      !  Apply bounds on ia
      if (ib .eq. 1) then
	iamin = 1
	iamax = N
      else
	iamin = ib + 1
	iamax = (N * ib) / (ib - 1)
	if (((iamax - ib) ** 2) .ge. N ) iamax = iamax + 1
      end if

      print *, "Running ib = ", ib, "ia from ", iamin, "to", iamax

      do ia = iamin, iamax

        !
        !  Check for common divisors and 
        !  only process box dimensions with none (above 1)
        !

        if (gcd(ia, ib) .eq. 1) then

	  box_area = ia * ib
	  !  max_area is the largest area for N boxes satisfying the area bound
	  max_area = box_area * (N + 1) - 1

	  min_area = box_area * ib

	  ! Bounds on ix and iy
	  ixmax = (N+1) * ib - 1
	  ixmin2 = ib * ib * ia
	  iymax2 = min(ixmax * ixmax, max_area)
	  iymin = ia

	  !
	  !  Calculate the best bounds on a/b
	  !  for a shrunken box

          call closest_boxes(ia, ib)

	  !
	  !  Generate a vector of possible pallet lengths
	  !

 	  call pallet_lengths(ia, ib, N)

	  !
 	  !  Initialise the msi tracking vectors
	  !
	  is_tight = 0
	  is_slack = 0
	  is_tight_low = 0
	  is_tight_high = 0

	  !
	  !  Iterate through all X and Y combinations
	  !
	  ix = izl(N_izl)
	  ix2 = ix * ix
	  iy = izl(1)
	  ix_index = N_izl
	  iy_index = 1

	  do while (ix2 .ge. ixmin2)

	    ipallet = ix * iy

	    if (ipallet .le. max_area .and. ipallet .ge. min_area) then

	      !
	      !  Check if we have a valid MSI and output if true
	      !
	      if (is_msi(ix, iy, ia, ib)) call output_msi(ix, iy, ia, ib)

	    end if

	    !
	    !  Move on to next combination of ix and iy
	    !
	    call increment_pallet(ix, iy, ix_index, iy_index)

	  end do

        end if

      end do

    end do

  end subroutine generate_all_msi



  !
  !  subroutine closest_boxes
  !
  !  This routine calculates the ratio of the closest smaller a/b pair
  !  relative to the input a/b box ratio
  !

  subroutine closest_boxes(ia, ib)

    implicit none
    integer(8), intent(in) :: ia, ib
    integer(8) i, j


    if (ib .eq. 1) then
      ial = ia
      iau = ia
      ibl = ib
      ibu = ib
      return
    end if

    if (ia .eq. ib) then
      ial = 1
      iau = 1
      ibl = 1
      ibu = 1
      return
    end if

    ial = 0
    ibl = 1
    iau = 1
    ibu = 0

    do i = (ib - 1), 1, -1
      ! integer division (floor)
      j = (i * ia) / ib

      if (j * ibl .gt. i * ial) then
 	ibl = i
	ial = j
      end if

      j = j + 1
      if (j * ibu .lt. i * iau) then
	ibu = i
	iau = j
      end if
    end do

  end subroutine closest_boxes



  !
  !  subroutine pallet_lengths
  !
  !  This routine generates the efficient set of possible pallet lengths
  !  with upto N boxes
  !

   subroutine pallet_lengths(ia, ib, N)

    implicit none
    integer(8), intent(in) :: ia, ib, N
    integer(8), dimension(1001000) :: ipace
    integer(8) i, j, k, ilength, max_len, idx, pt

    ! integer division
    max_len = 1 + (ixmax - 1) / 64
    if (max_len .gt. X_MAX2) then
      print *, "X_MAX2 too small. Recompile with larger value. Also check X_MAX", max_len, X_MAX2
      stop
    end if
    ipace = 0

    do i = 0, N
	k = min(ib, N-i)
	do j = 0, k

	    ilength = j * ia + i * ib
	    if (ilength .ge. iymin .and. ilength .le. ixmax) then
	        ! integer division
		idx = 1 +  (ilength - 1) / 64
		pt = mod(ilength, 64)
		ipace(idx) = ibset(ipace(idx), pt)
	    end if

	end do
    end do


    N_izl = 0
    izl = 0
    do i = 1, max_len
      do j = 1, 64

	  pt = ibits(ipace(i), j, 1)

	  if (pt .ne. 0) then
	    N_izl = N_izl + 1
	    if (N_izl .gt. X_MAX) then
	      print *, "X_MAX too small. Recompile with larger value"
	      stop
	    end if

	    izl(N_izl) = 64 * (i-1) + j
	  end if
      end do
    end do
 
  end subroutine pallet_lengths




  !
  !  function is_msi
  !
  !  Function checks if problem is a valid MSI and returns true if it is.
  !
  !  Uses a cached record of previous calculations
  !

  function is_msi(ix, iy, ia, ib)

    implicit none
    logical :: is_msi
    integer(8), intent(in) :: ix, iy, ia, ib
    integer(8) i, j
    logical, dimension(2) :: is_x, is_y

    is_msi = .TRUE.

    !
    ! if ib = 1, we have an MSI
    !
    if (ib .eq. 1) return

    !
    ! Start by checking if we have already found 
    ! iy or ix can't be reduced
    !
    
    ! integer division
    i = 1 + (iy - 1) / 64
    j = mod(iy, 64)
    if (ibits(is_tight(i), j, 1) .eq. 1) return

    if (ix .ne. iy) then

      ! first check if both bounds are violated 
      ! integer division
      i = 1 + (ix - 1) / 64
      j = mod(ix, 64)
      if (ibits(is_tight(i), j, 1) .eq. 1) return

    end if

    !
    ! Now we proceed to see if ix and iy in combination
    ! can be reduced. In both cases, we check if we've
    ! already processed them individually.
    ! is_slack indicates whether we have processed it
    ! is_tight_high/low indicates if either bound
    ! is violated.
    is_msi = .FALSE.
    is_x(1:2) = .FALSE.
    is_y(1:2) = .FALSE.

    !
    ! now evaluate ix
    !
    if (ibits(is_slack(i), j, 1) .ne. 0) then

      if (ibits(is_tight_low(i), j, 1)  .eq. 1) is_x(1) = .TRUE.
      if (ibits(is_tight_high(i), j, 1) .eq. 1) is_x(2) = .TRUE.

    else

      call check_shrinkable(ix, ia, ib, is_x)

      if (is_x(1) .and. is_x(2)) then

	is_tight(i) = ibset(is_tight(i), j)
        is_msi = .TRUE.
        return

      else

        if (is_x(1)) is_tight_low(i)  = ibset(is_tight_low(i), j)
        if (is_x(2)) is_tight_high(i) = ibset(is_tight_high(i), j)
	is_slack(i) = ibset(is_slack(i), j)

      end if

    end if 


    !
    ! If we have a square pallet we're done
    !
    if (ix .eq. iy) return


    !
    ! now evaluate iy
    !
    ! integer division
    i = 1 + (iy - 1) / 64
    j = mod(iy, 64)

    if (ibits(is_slack(i), j, 1) .ne. 0) then

      if (ibits(is_tight_low(i), j, 1)  .eq. 1) is_y(1) = .TRUE.
      if (ibits(is_tight_high(i), j, 1) .eq. 1) is_y(2) = .TRUE.

    else

      call check_shrinkable(iy, ia, ib, is_y)


      if (is_y(1) .and. is_y(2)) then

	is_tight(i) = ibset(is_tight(i), j)
	is_msi = .TRUE.
        return

      else

        if (is_y(1)) is_tight_low(i)  = ibset(is_tight_low(i), j)
        if (is_y(2)) is_tight_high(i) = ibset(is_tight_high(i), j)
	is_slack(i) = ibset(is_slack(i), j)
	is_msi = .FALSE.

      end if
 
    end if 

    !
    !  Finally, check if the opposite bounds are violated
    !
    if (is_x(1) .and. is_y(2)) is_msi = .TRUE.
    if (is_y(1) .and. is_x(2)) is_msi = .TRUE.


  end function is_msi




  !
  !  subroutine check_shrinkable
  !
  !  This checks if a pallet length iz can be shrunk when the
  !  dissections are determined by ia and ib
  !  returns the vector is_violated with .T. when the bound violates
  !  the dissections
  !

  subroutine check_shrinkable(iz, ia, ib, is_violated)

    implicit none
    integer(8), intent(in) :: iz, ia, ib
    logical, intent(inout) :: is_violated(2)
    integer(8) az, i, j
    integer(8) llb, lub, llbnew, ulb, uub, ulbnew, ibltemp, ibutemp

    az = iz / ia
    is_violated(1:2) = .FALSE.
    ibltemp = ibl * (ib - 1)
    ibutemp = ibu * (ib - 1)

    ! run through all the constraints to identify if any are
    ! violated for either upper or lower bounds.
    ! We can bail early if both are violated  
    do i = 0, az
      ! integer division...
      j = (iz - i * ia) / ib

      ! lower bound
      if (.not. is_violated(1)) then
	llbnew = ib * (i * ial + ibl * j)
	if (i .gt. 0) then
	  llb = max(llb, llbnew)
	  lub = min(lub, llbnew + ibltemp)
	else
	  llb = llbnew
	  lub = llbnew + ibltemp
	end if
	if (llb .gt. lub) is_violated(1) = .TRUE.
      end if

      ! upper bound
      if (.not. is_violated(2)) then
	ulbnew = ib * (i * iau + ibu * j)
	if (i .gt. 0) then
	  ulb = max(ulb, ulbnew)
	  uub = min(uub, ulbnew + ibutemp)
	else
	  ulb = ulbnew
	  uub = ulbnew + ibutemp
	end if

	if (ulb .gt. uub) is_violated(2) = .TRUE.	
      end if

      if (is_violated(1) .and. is_violated(2)) return

    end do

    ! check final constraint
    lub = min(lub, ((az + 1) * ib * ial - ibl))
    if (llb .gt. lub) is_violated(1) = .TRUE.

    uub = min(uub, ((az + 1) * ib * iau - ibu))
    if (ulb .gt. uub) is_violated(2) = .TRUE.	

  end subroutine check_shrinkable



  !
  !  subroutine increment_pallet
  !
  !  This handles incrementing the pallet dimensions within range
  !  it relies on the calling routine checking for ixmin .gt. ix as an indicator
  !  that we are done incrementing and its time to move on.
  !  start with big ix and small iy, increment iy and then decrement ix
  !

  subroutine increment_pallet(ix, iy, ix_index, iy_index)

    implicit none
    integer(8), intent(inout) :: ix, iy, ix_index, iy_index

    do while ((iy * iy) .le. iymax2)

      iy_index = iy_index + 1
      if (iy_index .le. N_izl) iy = izl(iy_index)

      ! check if we need to decrement ix
      if ((iy_index .gt. N_izl) .or. (iy .gt. ix) .or. ((iy * iy) .gt. iymax2) .or. ((ix * iy) .gt. max_area)) then

	! most likely this loop will only execute once...
inner:	do while (ix_index .gt. 0)

	  ix_index = ix_index - 1

	  if (ix_index .gt. 0) then

	    ix = izl(ix_index)
	    ix2 = ix * ix
	    iy_index = 1
	    iy = izl(iy_index)
	    
	    if (ix2 .lt. ixmin2) return

	    ! check that the max_area bound isn't violated
	    if ((ix * iy) .le. max_area) exit inner
	  else
	    ix2 = ixmin2 - 1
	    return
	  endif
	end do inner

      end if

      ! check that iy doesn't violate minimum bound before returning
      if (iy .ge. iymin) return

    end do

  end subroutine increment_pallet



