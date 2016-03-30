! ******************************************************************************
!
! $Id: readGridSCREAM.f90,v 1.2 2005/03/09 03:55:30 haselbac Exp $
!
! Filename: readGridSCREAM.F90
!
! Purpose: Read 2d grid file in SCREAM format.
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
! Copyright: (c) 2002 by the University of Illinois
!
! RCS Revision history:
!
!   $Log: readGridSCREAM.f90,v $
!   Revision 1.2  2005/03/09 03:55:30  haselbac
!   Changed formatted write statements
!
!   Revision 1.1  2005/03/05 17:25:53  haselbac
!   Initial revision
!
!   Revision 1.1  2003/03/07 15:26:34  haselbac
!   Initial revision
!
!
! ******************************************************************************

SUBROUTINE readGridSCREAM

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

  INTEGER :: bType,ib,ic,ie,iFile,iv
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: c2v
  CHARACTER :: choice
  CHARACTER*(MAX_STRING_LEN) :: dummyString,iFileName

! ******************************************************************************
! Start
! ******************************************************************************

  WRITE(STDOUT,'(/,1X,A)') 'Enter file name:'
  READ(STDIN,'(A)') iFileName
  WRITE(STDOUT,'(/)')

  iFile = FILE_UNIT_GRID_INPUT

  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)   
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(FILE_OPEN_ERROR,iFileName,errorFlag)
  END IF ! errorFlag

  WRITE(STDOUT,'(1X,A)') 'Reading grid file in SCREAM format...'
  
! ******************************************************************************
! Read grid file
! ******************************************************************************

  READ(iFile,*) grid%nCells,grid%nVert,grid%nBounds

  WRITE(STDOUT,'(3X,A,1X,I7)') 'Number of vertices:   ',grid%nVert
  WRITE(STDOUT,'(3X,A,1X,I7)') 'Number of cells:      ',grid%nCells
  WRITE(STDOUT,'(3X,A,1X,I7)') 'Number of boundaries: ',grid%nBounds  

! ==============================================================================
! Read coordinates 
! ==============================================================================

  ALLOCATE(grid%xy(2,grid%nVert),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%xy',errorFlag)
  END IF ! errorFlag
  
  WRITE(STDOUT,'(1X,A)') 'Coordinates...'
   
  DO iv = 1,grid%nVert 
    READ(iFile,*) grid%xy(1,iv),grid%xy(2,iv)
  END DO ! iv

! ==============================================================================
! Read connectivity
! ==============================================================================

  WRITE(STDOUT,'(1X,A)') 'Connectivity...'
 
  ALLOCATE(c2v(4,grid%nCells),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'c2v',errorFlag)
  END IF ! errorFlag

  DO ic = 1,grid%nCells
    READ(iFile,*) c2v(1,ic),c2v(2,ic),c2v(3,ic),c2v(4,ic)
  END DO ! ic

! ------------------------------------------------------------------------------
! Sort into triangles and quadrilaterals
! ------------------------------------------------------------------------------

  grid%nTris  = 0
  grid%nQuads = 0

  DO ic = 1,grid%nCells
    IF ( c2v(1,ic) == c2v(4,ic) ) THEN 
      grid%nTris = grid%nTris + 1
    ELSE 
      grid%nQuads = grid%nQuads + 1
    END IF ! c2v
  END DO ! ic

  WRITE(STDOUT,'(3X,A,1X,I7)') 'Number of triangles:      ',grid%nTris
  WRITE(STDOUT,'(3X,A,1X,I7)') 'Number of quadrilaterals: ',grid%nQuads
  
  IF ( grid%nTris > 0 ) THEN 
    ALLOCATE(grid%tri2v(3,grid%nTris),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'grid%tri2v',errorFlag)
    END IF ! errorFlag
  END IF ! grid%nTris

  IF ( grid%nQuads > 0 ) THEN 
    ALLOCATE(grid%quad2v(4,grid%nQuads),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'grid%quad2v',errorFlag)
    END IF ! errorFlag
  END IF ! grid%nQuads

  grid%nTris  = 0
  grid%nQuads = 0

  DO ic = 1,grid%nCells
    IF ( c2v(1,ic) == c2v(4,ic) ) THEN 
      grid%nTris = grid%nTris + 1
      grid%tri2v(1:3,grid%nTris) = c2v(1:3,ic)
    ELSE 
      grid%nQuads = grid%nQuads + 1
      grid%quad2v(1:4,grid%nQuads) = c2v(1:4,ic)      
    END IF ! c2v
  END DO ! ic

  DEALLOCATE(c2v,STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(DEALLOCATE_ERROR,'c2v',errorFlag)
  END IF ! errorFlag

! ==============================================================================
! Read boundary information
! ==============================================================================
 
  WRITE(STDOUT,'(1X,A)') 'Boundaries...' 
  
  ALLOCATE(grid%bound(grid%nBounds),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound',errorFlag)
  END IF ! errorFlag  
  
  DO ib = 1,grid%nBounds
    READ(iFile,*) bType,grid%bound(ib)%nEdges
    
    WRITE(STDOUT,'(3X,A,1X,I2)') 'Boundary:',ib
    WRITE(STDOUT,'(5X,A,1X,I2)') 'Type:',bType
    WRITE(STDOUT,'(5X,A,1X,I4)') 'Number of edges:',grid%bound(ib)%nEdges
    
    IF ( bType == 10 ) THEN ! SCREAM slip wall 
      grid%bound(ib)%bType = 400
      grid%bound(ib)%bName = 'Slip wall'
    ELSE IF ( bType == 11 ) THEN ! SCREAM no-slip wall (adiabatic)
      grid%bound(ib)%bType = 300
      grid%bound(ib)%bName = 'No-slip wall'      
    ELSE IF ( bType == 13 ) THEN ! SCREAM no-slip wall (isothermal) 
      grid%bound(ib)%bType = 310    
      grid%bound(ib)%bName = 'No-slip wall'      
    ELSE IF ( bType == 20 ) THEN ! SCREAM inflow
      grid%bound(ib)%bType = 100
      grid%bound(ib)%bName = 'Inflow'      
    ELSE IF ( bType == 30 ) THEN ! SCREAM outflow
      grid%bound(ib)%bType = 200
      grid%bound(ib)%bName = 'Outflow'      
    ELSE IF ( bType == 31 ) THEN ! SCREAM outflow
      grid%bound(ib)%bType = 200     
      grid%bound(ib)%bName = 'Outflow'        
    ELSE IF ( bType == 40 ) THEN ! SCREAM farfield
      grid%bound(ib)%bType = 500     
      grid%bound(ib)%bName = 'Farfield'         
    ELSE IF ( bType == 50 ) THEN ! SCREAM symmetry
      grid%bound(ib)%bType = 400 ! map to slip wall for now
      grid%bound(ib)%bName = 'Symmetry'      
    ELSE IF ( bType == 60 ) THEN ! SCREAM periodic
      grid%bound(ib)%bType = 400 ! map to slip wall for now
      grid%bound(ib)%bName = 'Periodic'      
    ELSE 
      CALL errorHandling(BTYPE_CONVERT_ERROR)                      
    END IF ! bType
    
    WRITE(STDOUT,'(5X,A,1X,I3)') 'New type:',grid%bound(ib)%bType
    
    ALLOCATE(grid%bound(ib)%e2v(2,grid%bound(ib)%nEdges),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'grid%bound%e2v',errorFlag)
    END IF ! errorFlag    
    
    DO ie = 1,grid%bound(ib)%nEdges
      READ(iFile,*) grid%bound(ib)%e2v(1,ie),grid%bound(ib)%e2v(2,ie)
    END DO ! ie
  END DO ! ib

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(iFile,IOSTAT=errorFlag)   
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(FILE_CLOSE_ERROR,iFileName,errorFlag)
  END IF ! errorFlag
   
  WRITE(STDOUT,'(1X,A,/)') 'Grid file read successfully.'
  
! ******************************************************************************
! End
! ******************************************************************************
  
END SUBROUTINE readGridSCREAM
