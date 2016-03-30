! ******************************************************************************
!
! Procedure: BinarySearchInteger.F90
!
! Purpose: Search array with elements in ascending order using binary search.
!
! Description: See appropriate textbooks.
!
! Input:
!   a   array with elements in ascending order
!   n   size of array
!   v   value which is being searched for
! 
! Output:
!   i   location of v in a
!
! Notes: 
!   1. If v is not in a, the error flag ELEMENT_NOT_FOUND is returned.
!
! Copyright: (c) 2000, 2001 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE BinarySearchInteger(a,n,v,i)

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
! Parameters
! ==============================================================================   
  
  INTEGER, INTENT(IN) :: n,v
  INTEGER, INTENT(IN) :: a(n)
  INTEGER, INTENT(OUT) :: i
  
! ==============================================================================
! Local variables
! ==============================================================================

  INTEGER :: il,im,iu
  
  
! ******************************************************************************
! Start
! ******************************************************************************

  il = 1 ! initialise lower limit
  iu = n ! initialise upper limit

  DO 
    im = (iu + il)/2 ! compute midpoint
    
    IF ( v < a(im) ) THEN 
      iu = im - 1 ! replace lower limit
    ELSE 
      il = im + 1 ! replace upper limit
    END IF ! v

    IF ( v == a(im) ) THEN ! element found
      i = im
      EXIT
    END IF ! v

    IF ( iu < il ) THEN ! element not found
      i = -1
      EXIT
    END IF ! iu 
  END DO ! <empty> 

! ******************************************************************************
! End
! ******************************************************************************
  
END SUBROUTINE BinarySearchInteger
