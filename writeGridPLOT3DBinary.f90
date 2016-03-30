! ******************************************************************************
!
! $Id: writeGridPLOT3DBinary.f90,v 1.2 2016/03/04 17:02:00 neal Exp $
!
! Filename: writeGridPLOT3DBinary.F90
!
! Purpose: Write 2d grid file in PLOT3D format.
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
! $Log: writeGridPLOT3DBinary.f90,v $
! Revision 1.2  2016/03/04 17:02:00  neal
! *** empty log message ***
!
! Revision 1.1  2016/03/04 16:46:55  neal
! Initial revision
!
! Revision 1.2  2005/02/18 20:38:46  haselbac
! Added reading of proper PLOT3D grids
!
! Revision 1.1  2004/12/27 15:35:05  haselbac
! Initial revision
!
! Revision 1.1  2004/12/27 15:05:24  haselbac
! Initial revision
!
! ******************************************************************************

SUBROUTINE writeGridPLOT3DBinary

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

  INTEGER :: dummyInteger,i,iFile,ix,iy,iz,nBlocks,nx,ny,nz
  CHARACTER :: choice
  CHARACTER*(MAX_STRING_LEN) :: iFileName

! ******************************************************************************
! Start
! ******************************************************************************

  WRITE(STDOUT,'(/,1X,A)') 'Enter file name:'
  READ(STDIN,'(A)') iFileName
  WRITE(STDOUT,'(/)')

  iFile = FILE_UNIT_GRID_INPUT

  OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="NEW",IOSTAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(FILE_OPEN_ERROR,iFileName,errorFlag)
  END IF ! errorFlag

  WRITE(STDOUT,'(1X,A)') 'Writing binary grid file in PLOT3D format...'

! ******************************************************************************
! Write grid file
! ******************************************************************************

  nBlocks = 1
  WRITE(iFile) nBlocks

  WRITE(STDOUT,'(3X,A,1X,I3)') 'nBlocks:',nBlocks 

  gridPLOT3D%nx = 10
  gridPLOT3D%ny = 20
  gridPLOT3D%nz = 30
  WRITE(iFile) gridPLOT3D%nx,gridPLOT3D%ny,gridPLOT3D%nz
  
  WRITE(STDOUT,'(3X,A,3(1X,I3))') 'Dimensions:',gridPLOT3D%nx, & 
                                                gridPLOT3D%ny, &
                                                gridPLOT3D%nz
  
  IF ( gridPLOT3D%nx == 1 .OR. gridPLOT3D%ny == 1 .OR. gridPLOT3D%nz == 1 ) THEN  
    gridPLOT3D%nCells = MAX(gridPLOT3D%nx-1,1)* & 
                        MAX(gridPLOT3D%ny-1,1)* & 
                        MAX(gridPLOT3D%nz-1,1)

    is2d = .TRUE.
  ELSE
    gridPLOT3D%nCells = (gridPLOT3D%nx-1)*(gridPLOT3D%ny-1)*(gridPLOT3D%nz-1)

    is2d = .FALSE.
  END IF ! gridPLOT3D%nx
  
  gridPLOT3D%nVert  = gridPLOT3D%nx*gridPLOT3D%ny*gridPLOT3D%nz

  WRITE(STDOUT,'(3X,A,1X,I6)') 'Number of vertices:',gridPLOT3D%nVert
  WRITE(STDOUT,'(3X,A,1X,I6)') 'Number of cells:',gridPLOT3D%nCells

! ==============================================================================
! Read coordinates
! ==============================================================================

  ALLOCATE(gridPLOT3D%xyz(3,gridPLOT3D%nx,gridPLOT3D%ny,gridPLOT3D%nz), & 
                          STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'gridPLOT3D%xyz',errorFlag)
  END IF ! errorFlag

  WRITE(STDOUT,'(3X,A)') 'Coordinates...'

! For AGARD 445.6 wing file from Lee-Rausch
!  DO i = 1,3
!    DO iz = 1,gridPLOT3D%nz
!      DO iy = 1,gridPLOT3D%ny
!        DO ix = 1,gridPLOT3D%nx
!          READ(iFile,*) gridPLOT3D%xyz(i,ix,iy,iz)                   
!        END DO ! ix
!      END DO ! iy
!    END DO ! iz
!  END DO ! i

! Initializing default coordinates
  DO iz = 1,gridPLOT3D%nz
    DO iy = 1,gridPLOT3D%ny
      DO ix = 1,gridPLOT3D%nx
        gridPLOT3D%xyz(1,ix,iy,iz) = ix                   
        gridPLOT3D%xyz(2,ix,iy,iz) = iy                   
        gridPLOT3D%xyz(3,ix,iy,iz) = iz                   
      END DO ! ix
    END DO ! iy
  END DO ! iz

  WRITE(iFile) ((((gridPLOT3D%xyz(i,ix,iy,iz),ix=1,gridPLOT3D%nx), & 
                                               iy=1,gridPLOT3D%ny), &   
                                               iz=1,gridPLOT3D%nz), & 
                                               i=1,3) 

! ******************************************************************************
! End
! ******************************************************************************

  CLOSE(iFile,IOSTAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(FILE_CLOSE_ERROR,iFileName,errorFlag)
  END IF ! errorFlag

  WRITE(STDOUT,'(1X,A,/)') 'Grid file written successfully.'

! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE writeGridPLOT3DBinary
