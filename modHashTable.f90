! ******************************************************************************
!
! $Id: modHashTable.f90,v 1.1 2005/03/05 17:25:53 haselbac Exp $
!
! Filename: modHashTable.F90
!
! Purpose: Collection of hashing-related procedures
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
! $Log: modHashTable.f90,v $
! Revision 1.1  2005/03/05 17:25:53  haselbac
! Initial revision
!
! Revision 1.1  2004/12/27 15:36:12  haselbac
! Initial revision
!
! Revision 1.1  2004/12/27 15:06:02  haselbac
! Initial revision
!
!
! ******************************************************************************

MODULE modHashTable

  USE modError
  USE modGlobals

  IMPLICIT NONE

  SAVE

! ******************************************************************************
! Data
! ******************************************************************************

! ==============================================================================
! Public
! ==============================================================================
    
  INTEGER, PUBLIC :: hashTableCollisions,hashTableSize
  INTEGER, ALLOCATABLE, DIMENSION(:), PUBLIC :: hashTable

! ==============================================================================
! Private
! ==============================================================================

  INTEGER, PARAMETER, PRIVATE :: HASHTABLE_INIT = 0 ! must be <= 0 
  INTEGER, PARAMETER, PRIVATE :: NPRIMES = 48 

  INTEGER, PRIVATE :: primeNumbers(NPRIMES) = &
             (/         1,        251,        379,        509, &
                      761,       1021,       1531,       2039, &
                     3067,       4093,       6143,       8191, &
                    12289,      16381,      24571,      32749, &
                    49139,      65521,      98297,     131071, &
                   196613,     262139,     393209,     524287, &
                   786431,    1048573,    1572853,    2097143, &
                  3145721,    4194301,    6291449,    8388593, &
                 12582917,   16777213,   25165807,   33554393, &
                 50331599,   67108859,  100663261,  134217689, &
                201326549,  268435399,  402653171,  536870909, &
                805306349, 1073741789, 1610612711, 2147483647   /)  


  CONTAINS
  
! ******************************************************************************
! Procedures
! ******************************************************************************
  
  
  
  
! ******************************************************************************
!
! Purpose: Create hash table.
!
! Description: None.
!
! Input:
!   size	Size of hash table
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE CreateHashTable(size)

    IMPLICIT NONE    

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: size

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,ih

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL FindNearestPrime(2*size,hashTableSize)      

! ******************************************************************************
!   Allocate memory and initialize
! ******************************************************************************

    ALLOCATE(hashTable(hashTableSize),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'hashTable',errorFlag)
    END IF ! errorFlag        

    DO ih = 1,hashTableSize
      hashTable(ih) = HASHTABLE_INIT
    END DO ! ih
    
    hashTableCollisions = 0      

  END SUBROUTINE CreateHashTable  






! ******************************************************************************
!
! Purpose: Destroy hash table
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
   
  SUBROUTINE DestroyHashTable

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************  

    DEALLOCATE(hashTable,STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR) THEN 
      CALL errorHandling(DEALLOCATE_ERROR,'hashTable',errorFlag)
    END IF ! errorFlag        

! ******************************************************************************
!   End
! ******************************************************************************  

  END SUBROUTINE DestroyHashTable






! ******************************************************************************
!
! Purpose: Find nearest prime for hash table size.
!
! Description: None.
!
! Input:
!   size	Original size
!
! Output: 
!   primeSize	Nearest prime to original size
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE FindNearestPrime(size,primeSize)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
!   Arguments
! ==============================================================================
      
    INTEGER, INTENT(IN) :: size
    INTEGER, INTENT(OUT) :: primeSize

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: i 

! ******************************************************************************
!   Start, find nearest prime
! ******************************************************************************

    DO i = 1,NPRIMES-1      
      IF ( size >= primeNumbers(i) .AND. size < primeNumbers(i+1) ) THEN 
        primeSize = primeNumbers(i+1)
        EXIT
      END IF ! size
    END DO ! i   

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE FindNearestPrime




! ******************************************************************************
!
! Purpose: Compute key from series of integers.
!
! Description: None.
!
! Input:
!   a		Set of integers
!   aSize	Size of set of integers
!
! Output: 
!   key		Key
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE HashBuildKey(a,aSize,key)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: aSize
    INTEGER, INTENT(IN) :: a(1:aSize)
    INTEGER, INTENT(OUT) :: key

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: i,term
    INTEGER, PARAMETER :: A_RAND = 31415, B_RAND = 27183

! ******************************************************************************
!   Start, compute key
! ******************************************************************************  


!    key = MOD(MOD(A_RAND*B_RAND,hashTableSize)*v1 + v2,HUGE(1)/10)
!    key = MOD(A_RAND*v1 + B_RAND*v2,hashTableSize) 

! - ABS function used to guard against negative keys from integer overflow   
!    key = ABS(ABS(A_RAND*v1) + ABS(MOD(A_RAND*B_RAND,hashTableSize)*v2))       

    term = A_RAND
    key  = 0

    DO i = 1,aSize      
      term = MOD(term*B_RAND,hashTableSize)

      IF ( term < 0 ) THEN 
        term = term + HUGE(1) + 1
      END IF ! term

      key  = MOD((term*key + a(i)),HUGE(1))

      IF ( key < 0 ) THEN 
        key = key + HUGE(1) + 1
      END IF ! key
    END DO ! i

! ******************************************************************************
!   End
! ******************************************************************************  

  END SUBROUTINE HashBuildKey






! ******************************************************************************
!
! Purpose: Hash face.
!
! Description: None.
!
! Input:
!   key         Key from which address is computed
!   icg         Global cell number
!   ifl         Local face number
!   fv          Face vertices (only three are needed)
!   nVert	Number of vertices in face (3 = triangle, 4 = quadrilateral)
!
! Output:
!   faceType    Flag indicating whether face is new or not
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE HashFace(iq,key,fv,nFaces,f2v,faceType)
  
    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iq,key
    INTEGER, INTENT(IN) :: fv(1:3)
    INTEGER, INTENT(OUT) :: faceType
    INTEGER, INTENT(INOUT) :: nFaces
    INTEGER, INTENT(INOUT) :: f2v(:,:)  

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: addr,collCntr,incr
  
! ******************************************************************************
!   Start
! ******************************************************************************
    
! ******************************************************************************
!   Construct address from key based on sorted face vertices
! ******************************************************************************
      
    CALL HashFuncPrimary(key,addr)  

! ******************************************************************************
!   Insert face into hash table
! ******************************************************************************  

    collCntr = 0

    DO    

! ==============================================================================
!     Entry not yet occupied
! ==============================================================================  

      IF ( hashTable(addr) == HASHTABLE_INIT ) THEN     
        faceType = 0

        nFaces = nFaces + 1     
        hashTable(addr) = nFaces

        f2v(1:3,nFaces) = fv(1:3)
        f2v(4  ,nFaces) = iq ! Mark face with old iq index

        EXIT

! ==============================================================================
!     Entry already occupied  
! ==============================================================================  

      ELSE 

! ------------------------------------------------------------------------------
!       Entry occupied by same face
! ------------------------------------------------------------------------------    

        IF ( fv(1) == f2v(1,hashTable(addr)) .AND. & 
             fv(2) == f2v(2,hashTable(addr)) .AND. & 
             fv(3) == f2v(3,hashTable(addr)) ) THEN
          faceType = 1 
          
          f2v(4,hashTable(addr)) = -1 ! Mark face as already existing  

          EXIT

! ------------------------------------------------------------------------------
!       Entry occupied by some other face: COLLISION
! ------------------------------------------------------------------------------    

        ELSE       
          hashTableCollisions = hashTableCollisions + 1

          IF ( collCntr == 0 ) THEN 
            CALL HashFuncSecondary(key,incr)
          END IF ! collCntr

          collCntr = collCntr + 1        
          addr = MOD(addr + incr,hashTableSize)

          IF ( addr == 0 ) THEN 
            addr = 1
          END IF ! addr   
        END IF ! fv

      END IF ! hashTable
    END DO ! <empty>
  
! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE HashFace




! ******************************************************************************
!
! Purpose: Primary hash function: Compute address from key value.
!
! Description: None.
!
! Input:
!   key         Key from which address is computed
!
! Output:
!   addr	Address
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE HashFuncPrimary(key,addr)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: key
    INTEGER, INTENT(OUT) :: addr

! ******************************************************************************
!   Start, compute address
! ******************************************************************************

!      addr = 1 ! simple test, will lead to many collisions, but works
    addr = 1 + MOD(key,hashTableSize)             

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE HashFuncPrimary







! ******************************************************************************
!
! Purpose: Secondary hash function: Compute address from key value
!
! Description: None.
!
! Input:
!   key         Key from which address is computed
!
! Output:
!   addr	Address
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE HashFuncSecondary(key,addr)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: key
    INTEGER, INTENT(OUT) :: addr

! ******************************************************************************
!   Start, compute address
! ******************************************************************************

!      addr = 1 ! simple test, will lead to many collisions, but works
    addr = 1 + MOD(key,hashTableSize-2)

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE HashFuncSecondary




END MODULE modHashTable
