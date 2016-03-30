! ******************************************************************************
!
! $Id: modSortSearch.f90,v 1.1 2005/03/05 17:25:53 haselbac Exp $
!
! Filename: modHashTable.F90
!
! Purpose: Collection of searching and sorting procedures
!
! Description:
!
! Input:
!
! Output:
!
! Notes:
!
! Author: Andreas Haselbacher
!
! Copyright: (c) 2004 by the University of Illinois
!
! RCS Revision history:
!
! $Log: modSortSearch.f90,v $
! Revision 1.1  2005/03/05 17:25:53  haselbac
! Initial revision
!
! Revision 1.1  2004/12/27 15:36:05  haselbac
! Initial revision
!
! Revision 1.1  2004/12/27 15:06:12  haselbac
! Initial revision
!
!
! ******************************************************************************

MODULE ModSortSearch

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
! Local variables
! ==============================================================================
 
  INTEGER, PARAMETER :: ELEMENT_NOT_FOUND = -1 ! must be <0 

! ******************************************************************************
! Module subroutines
! ******************************************************************************

  CONTAINS
  
  
  
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
      i = ELEMENT_NOT_FOUND
      EXIT
    END IF ! iu 
  END DO ! <empty> 

! ******************************************************************************
! End
! ******************************************************************************
  
END SUBROUTINE BinarySearchInteger
  
  
  
  
! ******************************************************************************
!
! Procedure: QuickSortInteger.F90
!
! Purpose: Sort elements of integer array a of size n into ascending order.
!
! Description: Uses quicksort method. 
!
! Input:
!   a   array containing unsorted elements
!   n   number of elements in array a
! 
! Output:
!   a   array containing sorted elements
!
! Notes:
!   1. Taken from WWW, cannot remember where... Seems to originate from 
!      Nicklaus Wirths book, see remark below.
!   2. No modifications to original, apart from a few cosmetic changes and 
!      from deletion of second array which was also being sorted along with a.
!
! Copyright: (c) 2000, 2001 by the University of Illinois
!
! ******************************************************************************

  SUBROUTINE QuickSortInteger(a,n)

! NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTHS PASCAL
! BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(INOUT) :: a(n)

! Local Variables

    INTEGER :: i, j, k, l, r, s, stackl(50), stackr(50)
    INTEGER :: w, x

    s = 1
    stackl(1) = 1
    stackr(1) = n

! KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.

 10 CONTINUE
    l = stackl(s)
    r = stackr(s)
    s = s - 1

! KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.

 20 CONTINUE
    i = l
    j = r
    k = (l+r) / 2
    x = a(k)

! REPEAT UNTIL I > J.

    DO
      DO
        IF (a(i).LT.x) THEN ! Search from lower end
          i = i + 1
          CYCLE
        ELSE
          EXIT
        END IF
      END DO

      DO
        IF (x.LT.a(j)) THEN ! Search from upper end
          j = j - 1
          CYCLE
        ELSE
          EXIT
        END IF
      END DO

      IF (i.LE.j) THEN ! Swap positions i & j
        w = a(i)
        a(i) = a(j)
        a(j) = w
        i = i + 1
        j = j - 1
        IF (i.GT.j) EXIT
      ELSE
        EXIT
      END IF
    END DO

    IF (j-l.GE.r-i) THEN
      IF (l.LT.j) THEN
        s = s + 1
        stackl(s) = l
        stackr(s) = j
      END IF
      l = i
    ELSE
      IF (i.LT.r) THEN
        s = s + 1
        stackl(s) = i
        stackr(s) = r
      END IF
      r = j
    END IF

    IF (l.LT.r) GO TO 20
    IF (s.NE.0) GO TO 10

    RETURN

  END SUBROUTINE QuickSortInteger





! ******************************************************************************
!
! Procedure: QuickSortIntegerInteger.F90
!
! Purpose: Sort elements of integer array a of size n into ascending order and 
!   sort b along with it.
!
! Description: Uses quicksort method. 
!
! Input:
!   a   array containing unsorted elements
!   b   integer array
!   n   number of elements in array a
! 
! Output:
!   a   array containing sorted elements
!   b   integer array
!
! Notes:
!   1. Taken from WWW, cannot remember where... Seems to originate from 
!      Nicklaus Wirths book, see remark below.
!   2. No modifications to original, apart from a few cosmetic changes and 
!      from deletion of second array which was also being sorted along with a.
!
! Copyright: (c) 2003 by the University of Illinois
!
! ******************************************************************************

  SUBROUTINE QuickSortIntegerInteger(a,b,n)

! NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTHS PASCAL
! BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(INOUT) :: a(n),b(n)

! Local Variables

    INTEGER :: i, j, k, l, r, s, stackl(50), stackr(50)
    INTEGER :: v, w, x

    s = 1
    stackl(1) = 1
    stackr(1) = n

! KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.

 10 CONTINUE
    l = stackl(s)
    r = stackr(s)
    s = s - 1

! KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.

 20 CONTINUE
    i = l
    j = r
    k = (l+r) / 2
    x = a(k)

! REPEAT UNTIL I > J.

    DO
      DO
        IF (a(i).LT.x) THEN ! Search from lower end
          i = i + 1
          CYCLE
        ELSE
          EXIT
        END IF
      END DO

      DO
        IF (x.LT.a(j)) THEN ! Search from upper end
          j = j - 1
          CYCLE
        ELSE
          EXIT
        END IF
      END DO
      
      IF (i.LE.j) THEN ! Swap positions i & j
        w = a(i)
        a(i) = a(j)
        a(j) = w
        v = b(i)
        b(i) = b(j)
        b(j) = v
        i = i + 1
        j = j - 1
        IF (i.GT.j) EXIT
      ELSE
        EXIT
      END IF
    END DO

    IF (j-l.GE.r-i) THEN
      IF (l.LT.j) THEN
        s = s + 1
        stackl(s) = l
        stackr(s) = j
      END IF
      l = i
    ELSE
      IF (i.LT.r) THEN
        s = s + 1
        stackl(s) = i
        stackr(s) = r
      END IF
      r = j
    END IF

    IF (l.LT.r) GO TO 20
    IF (s.NE.0) GO TO 10

    RETURN

  END SUBROUTINE QuickSortIntegerInteger


! ******************************************************************************
! End
! ******************************************************************************

END MODULE ModSortSearch

