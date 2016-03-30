! ******************************************************************************
!
! $Id: readGridPLOT3DBinary.f90,v 1.1 2016/03/04 16:46:09 neal Exp $
!
! Filename: readGridPLOT3DBinary.F90
!
! Purpose: Read 2d grid file in PLOT3D format.
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
! $Log: readGridPLOT3DBinary.f90,v $
! Revision 1.1  2016/03/04 16:46:09  neal
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

SUBROUTINE readGridPLOT3DBinary

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

  OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD",IOSTAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(FILE_OPEN_ERROR,iFileName,errorFlag)
  END IF ! errorFlag

  WRITE(STDOUT,'(1X,A)') 'Reading binary grid file in PLOT3D format...'

! ******************************************************************************
! Read grid file
! ******************************************************************************

! For AGARD 445.6 wing file from Lee-Rausch
!  READ(iFile,*) dummyInteger
!  READ(iFile,*) gridPLOT3D%nx
!  READ(iFile,*) gridPLOT3D%ny
!  READ(iFile,*) gridPLOT3D%nz

  READ(iFile) nBlocks

  WRITE(STDOUT,'(3X,A,1X,I3)') 'nBlocks:',nBlocks 

  READ(iFile) gridPLOT3D%nx,gridPLOT3D%ny,gridPLOT3D%nz
  
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

  READ(iFile) ((((gridPLOT3D%xyz(i,ix,iy,iz),ix=1,gridPLOT3D%nx), & 
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

  WRITE(STDOUT,'(1X,A,/)') 'Grid file read successfully.'

! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE readGridPLOT3DBinary
