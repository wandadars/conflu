! ******************************************************************************
!
! $Id: PLOT3D2Generic.f90,v 1.1 2005/03/05 17:25:53 haselbac Exp $
!
! Filename: PLOT3D2Generic.F90
!
! Purpose: Convert 2d grid in PLOT3D format into 2d generic grid.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes: 
!   1. Not general, does not work for grids with wake-cuts. 
!   2. Hard-coded for laminar flat plate with lead-in region.
!   3. Assumes that grid is 2d in x-y plane.
!
! Author: Andreas Haselbacher
!
! Copyright: (c) 2005 by the University of Illinois
!
! RCS Revision history:
!
! $Log: PLOT3D2Generic.f90,v $
! Revision 1.1  2005/03/05 17:25:53  haselbac
! Initial revision
!
! Revision 1.1  2005/02/18 20:37:52  haselbac
! Initial revision
!
! ******************************************************************************

SUBROUTINE PLOT3D2Generic

  USE modError
  USE modGlobals
  USE modGrid
  USE modSortSearch

  IMPLICIT NONE

  INTERFACE
    INTEGER FUNCTION ijk2l(i,j,k,ni,nj)
      INTEGER, INTENT(IN) :: i,j,k,ni,nj
    END FUNCTION
  END INTERFACE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Local variables
! ==============================================================================

  INTEGER :: i,ib,ie,iq,iv,j,k,ni,nj,offs

! ******************************************************************************
! Start, get additional information for extrusion
! ******************************************************************************

  WRITE(STDOUT,'(/,1X,A)') & 
    'Creating 2d-grid in generic format from PLOT3D grid...'

! ******************************************************************************
! Set up boundary data structure for 2d grid
! ******************************************************************************

! ==============================================================================
! Count number of boundaries
! ==============================================================================

  grid%nBounds = 5

  ALLOCATE(grid%bound(grid%nBounds),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound',errorFlag)
  END IF ! errorFlag

! ==============================================================================
! Specify number of edges and boundary type for 2d grid
! ==============================================================================

  ib = 1

  grid%bound(ib)%nEdges = (gridPLOT3D%nx-1)/4
  grid%bound(ib)%bName  = 'Symmetry'
  grid%bound(ib)%bType  = 500
  
  ib = 2

  grid%bound(ib)%nEdges = 3*(gridPLOT3D%nx-1)/4
  grid%bound(ib)%bName  = 'NoSlip'
  grid%bound(ib)%bType  = 300 
   
  ib = 3

  grid%bound(ib)%nEdges = (gridPLOT3D%ny-1)
  grid%bound(ib)%bName  = 'Outflow'
  grid%bound(ib)%bType  = 200 
  
  ib = 4

  grid%bound(ib)%nEdges = (gridPLOT3D%nx-1)
  grid%bound(ib)%bName  = 'Farfield'
  grid%bound(ib)%bType  = 800
  
  ib = 5

  grid%bound(ib)%nEdges = (gridPLOT3D%ny-1)
  grid%bound(ib)%bName  = 'Inflow'
  grid%bound(ib)%bType  = 100     

! ==============================================================================
! Specify edge connectivity for 2d grid
! ==============================================================================

  DO ib = 1,grid%nBounds
    ALLOCATE(grid%bound(ib)%e2v(2,grid%bound(ib)%nEdges),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN
      CALL errorHandling(ALLOCATE_ERROR,'grid%bound%e2v',errorFlag)
    END IF ! errorFlag
  END DO ! ib

  ni = grid%bound(1)%nEdges + grid%bound(2)%nEdges + 1
  nj = grid%bound(3)%nEdges + 1

  ib = 1
  
  DO ie = 1,grid%bound(ib)%nEdges
    grid%bound(ib)%e2v(1,ie) = ie
    grid%bound(ib)%e2v(2,ie) = ie + 1
  END DO ! ie
  
  ib = 2
  
  offs = grid%bound(1)%nEdges
  
  DO ie = 1,grid%bound(ib)%nEdges
    grid%bound(ib)%e2v(1,ie) = ie     + offs
    grid%bound(ib)%e2v(2,ie) = ie + 1 + offs
  END DO ! ie  

  ib = 3
      
  DO ie = 1,grid%bound(ib)%nEdges
    j = ie
  
    grid%bound(ib)%e2v(1,ie) = ijk2l(ni,j  ,1,ni,nj)    
    grid%bound(ib)%e2v(2,ie) = ijk2l(ni,j+1,1,ni,nj)
  END DO ! ie  

  ib = 4
            
  DO ie = grid%bound(ib)%nEdges,1,-1
    i = ie
  
    grid%bound(ib)%e2v(1,ie) = ijk2l(i+1,nj,1,ni,nj)    
    grid%bound(ib)%e2v(2,ie) = ijk2l(i  ,nj,1,ni,nj)
  END DO ! ie 
  
  ib = 5
      
  DO ie = grid%bound(ib)%nEdges,1,-1
    j = ie
  
    grid%bound(ib)%e2v(1,ie) = ijk2l(1,j+1,1,ni,nj)    
    grid%bound(ib)%e2v(2,ie) = ijk2l(1,j  ,1,ni,nj)
  END DO ! ie 
  
! ==============================================================================
! Store 2d connectivity
! ==============================================================================

  grid%nTris  = 0
  grid%nQuads = gridPLOT3D%nCells
  grid%nCells = grid%nTris + grid%nQuads

  ALLOCATE(grid%quad2v(4,grid%nQuads),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'grid%quad2v',errorFlag)
  END IF ! errorFlag

  ni = gridPLOT3D%nx
  nj = gridPLOT3D%ny  

  DO j = 1,nj-1
    DO i = 1,ni-1
      iq = ijk2l(i,j,1,ni-1,nj-1)
    
      grid%quad2v(1,iq) = ijk2l(i  ,j  ,1,ni,nj) 
      grid%quad2v(2,iq) = ijk2l(i+1,j  ,1,ni,nj)
      grid%quad2v(3,iq) = ijk2l(i+1,j+1,1,ni,nj)
      grid%quad2v(4,iq) = ijk2l(i  ,j+1,1,ni,nj)
    END DO ! i
  END DO ! j

! ==============================================================================
! Store 2d coordinates
! ==============================================================================

  grid%nVert = gridPLOT3D%nVert

  ALLOCATE(grid%xy(2,grid%nVert),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'grid%xy',errorFlag)
  END IF ! errorFlag

  ni = gridPLOT3D%nx
  nj = gridPLOT3D%ny  

  DO j = 1,nj
    DO i = 1,ni
      iv = ijk2l(i,j,1,ni,nj) 
    
      grid%xy(1,iv) = gridPLOT3D%xyz(1,i,j,1) 
      grid%xy(2,iv) = gridPLOT3D%xyz(2,i,j,1)
    END DO ! i
  END DO ! j

! ******************************************************************************
! Destroy PLOT3D data 
! ******************************************************************************

  DEALLOCATE(gridPLOT3D%xyz,STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(DEALLOCATE_ERROR,'gridPLOT3D%xyz',errorFlag)
  END IF ! errorFlag
  
! ******************************************************************************
! End
! ******************************************************************************

  WRITE(STDOUT,'(1X,A)') & 
    'Creating 2d-grid in generic format from PLOT3D grid done.'

END SUBROUTINE PLOT3D2Generic
