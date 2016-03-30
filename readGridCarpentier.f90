! ******************************************************************************
!
! $Id: readGridCarpentier.f90,v 1.1 2005/03/05 17:25:53 haselbac Exp $
!
! Filename: readGridCarpentier.F90
!
! Purpose: Read 2d grid file in Romuald Carpentiers format.
!
! Description: None.
!
! Input: None.
! 
! Output: None.
!
! Notes: 
!   1. This routine assumes that the grid is triangular.
!   2. This routine assumes that the grid has been obtained from a 
!      structurded quadrilateral grid by inserting diagonal edges.
!   3. The number of x- and y-grid points has been hard-coded!
!   4. For conversion to quadrilateral grid, assume that the 
!      triangular cells are always defined in the same way!
!
! Author: Andreas Haselbacher
!
! Copyright: (c) 2001 by the University of Illinois
!
! RCS Revision history:
!
!   $Log: readGridCarpentier.f90,v $
!   Revision 1.1  2005/03/05 17:25:53  haselbac
!   Initial revision
!
!   Revision 1.1  2003/03/07 15:26:34  haselbac
!   Initial revision
!
!
! ******************************************************************************

SUBROUTINE readGridCarpentier

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

  INTEGER :: i,ib,ie,iFile,iq,it,j,nx,ny
  CHARACTER :: choice
  CHARACTER*(MAX_STRING_LEN) :: iFileName

!  TEMPORARY for ROCFLO
!   INTEGER :: ix,iy,iz,nz
!   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: xr,yr,zr
!  END TEMPORARY

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

  WRITE(STDOUT,'(1X,A)') 'Reading grid file in Carpentier format...'
  
! ******************************************************************************
! Read grid file
! ******************************************************************************

  READ(iFile,*) grid%nVert,grid%nCells

  grid%nTris  = grid%nCells
  grid%nQuads = 0
  
  WRITE(STDOUT,'(3X,A,1X,I6)') 'Number of vertices:',grid%nVert
  WRITE(STDOUT,'(3X,A,1X,I6)') 'Number of triangles:',grid%nTris

  nx = 318 ! NOTE hard-code
  ny = 31 

! ==============================================================================
! Read coordinates 
! ==============================================================================

  ALLOCATE(grid%xy(2,grid%nVert),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%xyz',errorFlag)
  END IF ! errorFlag
  
  WRITE(STDOUT,'(3X,A)') 'Coordinates...'
   
  READ(iFile,*) ((grid%xy(i,j),i=1,2),j=1,grid%nVert)
  
! ==============================================================================
! Read connectivity
! ==============================================================================
 
  ALLOCATE(grid%tri2v(3,grid%nTris),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%tri2v',errorFlag)
  END IF ! errorFlag

  WRITE(STDOUT,'(3X,A)') 'Connectivity...'

  READ(iFile,*) ((grid%tri2v(i,j),i=1,3),j=1,grid%nTris)

! ==============================================================================
! Convert to quadrilateral mesh if desired
! ==============================================================================

  WRITE(STDOUT,'(3X,A)') 'Convert to quadrilateral grid? (y/n)'
  READ(STDIN,*) choice

  IF ( choice == 'y' ) THEN 
    grid%nQuads = grid%nTris/2
    grid%nCells = grid%nQuads

    ALLOCATE(grid%quad2v(4,grid%nQuads),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'grid%quad2v',errorFlag)
    END IF ! errorFlag

    iq = 0

    DO it = 1,grid%nTris,2
      iq = iq + 1
      grid%quad2v(1,iq) = grid%tri2v(3,it)
      grid%quad2v(2,iq) = grid%tri2v(1,it+1)
      grid%quad2v(3,iq) = grid%tri2v(1,it)
      grid%quad2v(4,iq) = grid%tri2v(2,it)
    END DO ! it    

    grid%nTris = 0
 
    DEALLOCATE(grid%tri2v,STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(DEALLOCATE_ERROR,'grid%tri2v',errorFlag)
    END IF ! errorFlag
    
    WRITE(STDOUT,'(3X,A,1X,I6)') 'Number of triangles:',grid%nTris
    WRITE(STDOUT,'(3X,A,1X,I6)') 'Number of quadrilaterals:',grid%nQuads
    WRITE(STDOUT,'(3X,A,1X,I6)') 'Number of cells:',grid%nCells    
  END IF ! choice

! ****************************************************************************** 
! Set remaining quantities
! ******************************************************************************

  grid%nBounds = 5

  ALLOCATE(grid%bound(grid%nBounds),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound',errorFlag)
  END IF ! errorFlag

! ------------------------------------------------------------------------------
! Boundary 1 - head end solid wall
! ------------------------------------------------------------------------------ 
 
  grid%bound(1)%bName = 'Head end wall'
  grid%bound(1)%bType = 300 
  grid%bound(1)%nEdges = ny-1 

  ALLOCATE(grid%bound(1)%e2v(2,grid%bound(1)%nEdges),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound%e2v',errorFlag)
  END IF ! errorFlag

  DO ie = 1,grid%bound(1)%nEdges
    grid%bound(1)%e2v(1,ie) = nx*(ny-1)+1 - (ie-1)*nx
    grid%bound(1)%e2v(2,ie) = grid%bound(1)%e2v(1,ie) - nx
  END DO ! ie

! ------------------------------------------------------------------------------
! Boundary 2 - injection surface - NOTE further hard-code
! ------------------------------------------------------------------------------ 

  grid%bound(2)%bName = 'Injection surface'
  grid%bound(2)%bType = 300 
  grid%bound(2)%nEdges = 83   

  ALLOCATE(grid%bound(2)%e2v(2,grid%bound(2)%nEdges),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound%e2v',errorFlag)
  END IF ! errorFlag

  DO ie = 1,grid%bound(2)%nEdges
    grid%bound(2)%e2v(1,ie) = ie
    grid%bound(2)%e2v(2,ie) = grid%bound(2)%e2v(1,ie) + 1  
  END DO ! ie

! ------------------------------------------------------------------------------
! Boundary 3 - solid wall - NOTE further hard code
! ------------------------------------------------------------------------------ 

  grid%bound(3)%bName = 'Casing wall'
  grid%bound(3)%bType = 300 
  grid%bound(3)%nEdges = nx-1 - grid%bound(2)%nEdges

  ALLOCATE(grid%bound(3)%e2v(2,grid%bound(3)%nEdges),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound%e2v',errorFlag)
  END IF ! errorFlag

  DO ie = 1,grid%bound(3)%nEdges
    grid%bound(3)%e2v(1,ie) = grid%bound(2)%nEdges + ie
    grid%bound(3)%e2v(2,ie) = grid%bound(3)%e2v(1,ie) + 1 
  END DO ! ie

! ------------------------------------------------------------------------------
! Boundary 4 - outlet 
! ------------------------------------------------------------------------------ 

  grid%bound(4)%bName = 'Outlet'
  grid%bound(4)%bType = 200 
  grid%bound(4)%nEdges = ny-1 

  ALLOCATE(grid%bound(4)%e2v(2,grid%bound(4)%nEdges),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound%e2v',errorFlag)
  END IF ! errorFlag

  DO ie = 1,grid%bound(4)%nEdges
    grid%bound(4)%e2v(1,ie) = ie*nx 
    grid%bound(4)%e2v(2,ie) = grid%bound(4)%e2v(1,ie) + nx
  END DO ! ie

! ------------------------------------------------------------------------------
! Boundary 5 - symmetry boundary
! ------------------------------------------------------------------------------ 

  grid%bound(5)%bName = 'Symmetry'
  grid%bound(5)%bType = 500 
  grid%bound(5)%nEdges = nx-1   

  ALLOCATE(grid%bound(5)%e2v(2,grid%bound(5)%nEdges),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound%e2v',errorFlag)
  END IF ! errorFlag

  DO ie = 1,grid%bound(5)%nEdges
    grid%bound(5)%e2v(1,ie) = nx*ny - (ie-1)
    grid%bound(5)%e2v(2,ie) = grid%bound(5)%e2v(1,ie) - 1 
  END DO ! ie

! BEGIN DEBUG
!  DO ib = 1,grid%nBounds
!    WRITE(STDOUT,*) '*** ',ib,grid%bound(ib)%nEdges
!    DO ie = 1,grid%bound(ib)%nEdges
!      WRITE(STDOUT,*) grid%bound(ib)%e2v(1,ie),grid%bound(ib)%e2v(2,ie)
!    END DO ! ie
!  END DO ! ib
! END DEBUG

! ******************************************************************************
! Write ROCFLO file - temporary
! ******************************************************************************
!
!  nz = 90
!
!  ALLOCATE(xr(nx,ny,nz),yr(nx,ny,nz),zr(nx,ny,nz),STAT=errorFlag)
!  IF ( errorFlag /= NO_ERROR ) THEN 
!    CALL errorHandling(ALLOCATE_ERROR,'xr,yr,zr',errorFlag)
!  END IF ! errorFlag
!
!  DO ix = 1,nx
!    DO iy = 1,ny
!      DO iz = 1,nz
!        xr(ix,iy,iz) = grid%xy(1,ix+(iy-1)*nx)
!        yr(ix,iy,iz) = grid%xy(2,ix+(iy-1)*nx)
!	zr(ix,iy,iz) = (iz-1)*0.141D0/(1.0D0*(nz-1))	
!      END DO ! iz
!    END DO ! iy
!  END DO ! ix
!
!  iFileName = 'onera_c1_001.msh'
!  OPEN(45,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN",IOSTAT=errorFlag)
!  IF ( errorFlag /= NO_ERROR ) THEN 
!    CALL errorHandling(FILE_OPEN_ERROR,iFileName,errorFlag)
!  END IF ! errorFlag
!
!  WRITE(45,*) 0.0
!  WRITE(45,*) nx,ny,nz
!  WRITE(45,*) (((xr(ix,iy,iz),yr(ix,iy,iz),zr(ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz)
!
!  CLOSE(45,IOSTAT=errorFlag)
!  IF ( errorFlag /= NO_ERROR ) THEN 
!    CALL errorHandling(FILE_CLOSE_ERROR,iFileName,errorFlag)
!  END IF ! errorFlag
!
!  DEALLOCATE(xr,yr,zr,STAT=errorFlag)
!  IF ( errorFlag /= NO_ERROR ) THEN 
!    CALL errorHandling(DEALLOCATE_ERROR,'xr,yr,zr',errorFlag)
!  END IF ! errorFlag

! ******************************************************************************
! End  
! ******************************************************************************

  CLOSE(iFile,IOSTAT=errorFlag)   
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(FILE_CLOSE_ERROR,iFileName,errorFlag)
  END IF ! errorFlag
   
  WRITE(STDOUT,'(1X,A,/)') 'Grid file read successfully.'
  
! ******************************************************************************
! End
! ******************************************************************************
  
END SUBROUTINE readGridCarpentier
