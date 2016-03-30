! ******************************************************************************
!
! $Id: distortVertices.f90,v 1.7 2016/03/04 16:42:23 neal Exp $
!
! Filename: distortVertices.F90
!
! Purpose: Distort interior vertices.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes: None.
!
! Author: Andreas Haselbacher
!
! Copyright: (c) 2004 by the University of Illinois
!
! RCS Revision history:
!
! $Log: distortVertices.f90,v $
! Revision 1.7  2016/03/04 16:42:23  neal
! Make error checking compliant with how the rest of the code handles errors.
!
! Revision 1.6  2016/03/04 16:35:53  neal
! *** empty log message ***
!
! Revision 1.5  2016/03/04 16:35:19  neal
! *** empty log message ***
!
! Revision 1.4  2016/03/04 16:34:15  neal
! *** empty log message ***
!
! Revision 1.3  2016/03/04 16:33:32  neal
! *** empty log message ***
!
! Revision 1.2  2016/03/04 16:31:31  neal
! Added a check for the required number of seeds. This is different for some compilers.
!
! Revision 1.1  2005/03/05 17:25:53  haselbac
! Initial revision
!
! Revision 1.2  2004/02/05 18:20:21  haselbac
! Determine grid spacing (on boundaries) automatically
!
! Revision 1.1  2004/02/05 17:43:34  haselbac
! Initial revision
!
! ******************************************************************************

SUBROUTINE distortVertices

  USE modError
  USE modGlobals
  USE modGrid

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Local variables
! ==============================================================================

  INTEGER :: flag,ib,ie,iv,v1,v2
  INTEGER :: randomSeed(2),nSeeds,I
  INTEGER, ALLOCATABLE :: TmpRandomSeeds(:)
  DOUBLE PRECISION :: distFrac,ds,dsMin,dx,dy,rnd
  CHARACTER(1) :: choice

! ******************************************************************************
! Start
! ******************************************************************************

  randomSeed(1) = 6971
  randomSeed(2) = 3029

  CALL RANDOM_SEED(SIZE = nSeeds)

  ALLOCATE ( TmpRandomSeeds(nSeeds), STAT=errorFlag)
  IF (errorFlag/= NO_ERROR ) THEN
     CALL errorHandling(ALLOCATE_ERROR,'TmpRandomSeeds',errorFlag)
  END IF   

  IF(nSeeds == 2) THEN

    TmpRandomSeeds(1) = randomSeed(1)
    TmpRandomSeeds(2) = randomSeed(2)

  ELSE
    ! We'll just populate this with odd seeds = Seed(1) and even seeds = Seed(2)
    DO I = 1,nSeeds,2
       TmpRandomSeeds(I) = randomSeed(1)
       END DO
  
    DO I = 2,nSeeds,2
       TmpRandomSeeds(I) = randomSeed(2)
    END DO
                     
                  
    CALL RANDOM_SEED ( PUT=TmpRandomSeeds )
    DEALLOCATE(TmpRandomSeeds, STAT=errorFlag)
    IF (errorFlag/= NO_ERROR ) THEN
       CALL errorHandling(DEALLOCATE_ERROR,'TmpRandomSeeds',errorFlag)
    END IF 
                           
  END IF

  WRITE(STDOUT,'(1X,A)') 'Do you want to distort grid? (y/n)'
  READ(STDIN,'(A)') choice
  
  IF ( choice .EQ. 'y' ) THEN 
    WRITE(STDOUT,'(1X,A)') 'Determining minimum spacing...'
  
    dsMin = HUGE(1.0D0)
  
    DO ib = 1,grid%nBounds
      DO ie = 1,grid%bound(ib)%nEdges
        v1 = grid%bound(ib)%e2v(1,ie)
        v2 = grid%bound(ib)%e2v(2,ie)

        dx = grid%xy(1,v2) - grid%xy(1,v1)
        dy = grid%xy(2,v2) - grid%xy(2,v1)
        ds = SQRT(dx*dx + dy*dy)

        dsMin = MIN(dsMin,ds)
      END DO ! ie
    END DO ! ib    
  
    WRITE(STDOUT,*) 'Minimum spacing:',dsMin
    WRITE(STDOUT,'(1X,A)') 'Enter distortion fraction:'
    READ(STDIN,*) distFrac
  
    DO iv = 1,grid%nVert
      flag = 0 

      boundLoop: DO ib = 1,grid%nBounds
        DO ie = 1,grid%bound(ib)%nEdges
          v1 = grid%bound(ib)%e2v(1,ie)
          v2 = grid%bound(ib)%e2v(2,ie)

          IF ( iv == v1 .OR. iv == v2 ) THEN 
            flag = 1
            EXIT boundLoop
          END IF ! iv
        END DO ! ie
      END DO boundLoop 

      IF ( flag == 0 ) THEN ! iv is not on boundaries 
        CALL RANDOM_NUMBER(rnd)      
        grid%xy(1,iv) = grid%xy(1,iv) + (2.0D0*rnd - 1.0D0)*distFrac*dsMin

        CALL RANDOM_NUMBER(rnd)
        grid%xy(2,iv) = grid%xy(2,iv) + (2.0D0*rnd - 1.0D0)*distFrac*dsMin          
      END IF ! flag
    END DO ! iv
  END IF ! choice

! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE distortVertices
