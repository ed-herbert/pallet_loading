
!
!   misc_functions.f95
!
!   Generates all minimum size instances of the pallet loading problem for
!   less than N+1 boxes (by area).
!
!   Miscellaneous functions
!
!   For details on the algorithm, see 
!   "An algorithm for identifying equivalence classes of the pallet loading problem"
!
!   Ed Herbert, 12 Nov 2021
!
!




  !
  !  function gcd
  !  Returns the greatest common divisor of integers i and j
  !

  function gcd(i, j)

    implicit none
    integer(8), intent(in)  :: i, j
    integer(8) gcd
    integer(8) d, t

    gcd = i
    d = j

    do while (d .ne. 0)

        t = d
        d = mod(gcd, d)
        gcd = t

    end do


  end function gcd

