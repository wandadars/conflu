! ******************************************************************************
!
! $Id: readGridCOBALT.f90,v 1.1 2005/03/05 17:25:53 haselbac Exp $
!
! Filename: readGridCOBALT.F90
!
! Purpose: Read grid in COBALT binary double precision format.
!
! Description: None.
!
! Input: None.
! 
! Output: None.
!
! Notes: 
!   1. Restricted to three-dimensional COBALT grids.
!   2. Restricted to single-zone COBALT grids.
!
! Author: Andreas Haselbacher
!
! Copyright: (c) 2001 by the University of Illinois
!
! RCS Revision history:
!
!   $Log: readGridCOBALT.f90,v $
!   Revision 1.1  2005/03/05 17:25:53  haselbac
!   Initial revision
!
!   Revision 1.1  2003/03/07 15:26:34  haselbac
!   Initial revision
!
!
! ******************************************************************************

SUBROUTINE readGridCOBALT

  USE modError
  USE modGlobals
  USE modGrid

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
! Parameters
! ==============================================================================   
  
! ==============================================================================
! Local variables
! ==============================================================================
 
  INTEGER :: i,iFile,j,nVertPerFace
  CHARACTER*(MAX_STRING_LEN) :: iFileName

! ******************************************************************************
! Start
! ******************************************************************************

  gridCENTAUR%nTets = 0
  gridCENTAUR%nHexs = 0
  gridCENTAUR%nPris = 0
  gridCENTAUR%nPyrs = 0

  gridCOBALT%nTris  = 0
  gridCOBALT%nQuads = 0
 
  WRITE(STDOUT,'(1X,A)') 'Enter file name:'
  READ(STDIN,'(A)') iFileName
  WRITE(STDOUT,'(/)')

  iFile     = FILE_UNIT_GRID_INPUT
 
  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", & 
       IOSTAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(FILE_OPEN_ERROR,iFileName,errorFlag)
  END IF ! errorFlag
  
! ******************************************************************************
! Read dimensions
! ******************************************************************************

  WRITE(STDOUT,'(/,1X,A)') 'Reading COBALT grid file...'
  WRITE(STDOUT,'(3X,A)') 'Dimensions...'

  READ(iFile,*) gridCOBALT%nDim,gridCOBALT%nZones,gridCENTAUR%nBounds

  IF ( gridCOBALT%nDim /= 3 ) THEN 
    CALL errorHandling(GRID_COBALT_NDIM_ERROR)
  END IF ! gridCOBALT

  IF ( gridCOBALT%nZones /= 1 ) THEN 
    CALL errorHandling(GRID_COBALT_NZONES_ERROR)
  END IF ! gridCOBALT

  READ(iFile,*) gridCENTAUR%nVert,gridCOBALT%nFaces,gridCENTAUR%nCells, & 
                gridCOBALT%nVertPerFaceMax,gridCOBALT%nFacesPerCellMax

! ******************************************************************************
! Allocate memory for coordinates and connectivity
! ******************************************************************************

  ALLOCATE(gridCENTAUR%xyz(3,gridCENTAUR%nVert),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%xyz')
  END IF ! errorFlag

  ALLOCATE(gridCOBALT%f2v(gridCOBALT%nVertPerFaceMax,gridCOBALT%nFaces), & 
           STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'gridCOBALT%f2v')
  END IF ! errorFlag

  gridCOBALT%f2v(:,:) = 0

  ALLOCATE(gridCOBALT%f2c(2,gridCOBALT%nFaces),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'gridCOBALT%f2c')
  END IF ! errorFlag

  ALLOCATE(gridCOBALT%nvpf(gridCOBALT%nFaces),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'gridCOBALT%nvpf')
  END IF ! errorFlag

! ******************************************************************************
! Read coordinates and connectivity 
! ******************************************************************************

  WRITE(STDOUT,'(3X,A)') 'Coordinates...'

  DO i = 1,gridCENTAUR%nVert
    READ(iFile,*) (gridCENTAUR%xyz(j,i),j=1,3) 
  END DO ! i

  WRITE(STDOUT,'(3X,A)') 'Connectivity...'

  DO i = 1,gridCOBALT%nFaces
    READ(iFile,*) nVertPerFace
    BACKSPACE(iFile)

    IF ( nVertPerFace == 3 ) THEN
      gridCOBALT%nTris = gridCOBALT%nTris + 1
    ELSE IF ( nVertPerFace == 4 ) THEN
      gridCOBALT%nQuads = gridCOBALT%nQuads + 1
    ELSE  
      CALL errorHandling(INVALID_FACETYPE_ERROR)
    END IF ! nVertPerFace

    READ(iFile,*) gridCOBALT%nvpf(i),(gridCOBALT%f2v(j,i),j=1,nVertPerFace), & 
                 (gridCOBALT%f2c(j,i),j=1,2)
  END DO ! i 

  CLOSE(iFile,IOSTAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(FILE_CLOSE_ERROR,iFileName,errorFlag)
  END IF ! errorFlag 

  WRITE(STDOUT,'(1X,A)') 'File read successfully.'

! ******************************************************************************
! Print grid statistics 
! ******************************************************************************

  WRITE(STDOUT,'(/,1X,A)') 'Grid statistics:'
  WRITE(STDOUT,'(3X,A,3X,I9)') 'Vertices:           ',gridCENTAUR%nVert
  WRITE(STDOUT,'(3X,A,3X,I9)') 'Cells:              ',gridCENTAUR%nCells
  WRITE(STDOUT,'(3X,A,3X,I9)') 'Faces:              ',gridCOBALT%nFaces
  WRITE(STDOUT,'(5X,A,1X,I9)') 'Triangular faces:   ',gridCOBALT%nTris
  WRITE(STDOUT,'(5X,A,1X,I9)') 'Quadrilateral faces:',gridCOBALT%nQuads
  WRITE(STDOUT,'(3X,A,3X,I9)') 'Boundary patches:   ',gridCENTAUR%nBounds


! ******************************************************************************
! End
! ******************************************************************************
  
END SUBROUTINE readGridCOBALT
