! ******************************************************************************
!
! $Id: readGridCENTAUR.f90,v 1.2 2005/03/07 00:56:30 haselbac Exp $
!
! Filename: readGridCENTAUR.F90
!
! Purpose: Read grid file from CENTAUR in binary double-precision format.
!
! Description:
!
! Input:
! 
! Output:
!
! Notes:
!   1. Cell and node pointers are written as 0.
!
! Author: Andreas Haselbacher
!
! Copyright: (c) 2004 by the University of Illinois
!
! RCS Revision history:
!
!   $Log: readGridCENTAUR.f90,v $
!   Revision 1.2  2005/03/07 00:56:30  haselbac
!   Added setting of nCells
!
!   Revision 1.1  2005/03/05 17:25:53  haselbac
!   Initial revision
!
!   Revision 1.2  2004/07/19 18:55:55  haselbac
!   Added writing of grid info
!
!   Revision 1.1  2004/07/17 22:14:21  haselbac
!   Initial revision
!
!
! ******************************************************************************

SUBROUTINE readGridCENTAUR

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

  INTEGER :: i,idum,iFile,j,quadOffs,triOffs
  CHARACTER(MAX_STRING_LEN) :: caseName,iFileName

! ******************************************************************************
! Open file and read title
! ******************************************************************************

  WRITE(STDOUT,'(/,1X,A)') 'Reading CENTAUR grid file...'
  WRITE(STDOUT,'(3X,A)') 'Enter case name:'
  READ(STDIN,*) caseName

  iFile     = FILE_UNIT_GRID_INPUT
  iFileName = TRIM(caseName)//'.hyb.bin'

  OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD",IOSTAT=errorFlag)   
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(FILE_OPEN_ERROR,iFileName,errorFlag)
  END IF ! errorFlag

  READ(iFile) gridCENTAUR%title
  
! ==============================================================================
! Coordinates
! ==============================================================================

  READ(iFile) gridCENTAUR%nVert
  
  WRITE(STDOUT,'(3X,A)') 'Coordinates...'
   
  ALLOCATE(gridCENTAUR%xyz(3,gridCENTAUR%nVert),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%xyz',errorFlag)
  END IF ! errorFlag 
   
  READ(iFile) ((gridCENTAUR%xyz(i,j),j=1,gridCENTAUR%nVert),i=1,3)  

  idum = 0
  READ(iFile) (idum,i=1,gridCENTAUR%nVert)

! ==============================================================================
! Cell connectivity
! ==============================================================================

  READ(iFile) gridCENTAUR%nTets
  IF ( gridCENTAUR%nTets > 0 ) THEN 
    ALLOCATE(gridCENTAUR%tet2v(4,gridCENTAUR%nTets),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%tet2v',errorFlag)
    END IF ! errorFlag    
  
    WRITE(STDOUT,'(3X,A)') 'Tetrahedra...'
    READ(iFile) ((gridCENTAUR%tet2v(i,j),j=1,gridCENTAUR%nTets),i=1,4)
    READ(iFile) (idum,i=1,gridCENTAUR%nTets)    
  END IF ! nTetrahedra
  
  READ(iFile) gridCENTAUR%nHexs
  IF ( gridCENTAUR%nHexs > 0 ) THEN 
    ALLOCATE(gridCENTAUR%hex2v(8,gridCENTAUR%nHexs),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%hex2v',errorFlag)
    END IF ! errorFlag  
  
    WRITE(STDOUT,'(3X,A)') 'Hexahedra...' 
    READ(iFile) ((gridCENTAUR%hex2v(i,j),j=1,gridCENTAUR%nHexs),i=1,8)
    READ(iFile) (idum,i=1,gridCENTAUR%nHexs)   
  END IF ! nHexahedra

  READ(iFile) gridCENTAUR%nPris
  IF ( gridCENTAUR%nPris > 0 ) THEN
    ALLOCATE(gridCENTAUR%pri2v(6,gridCENTAUR%nPris),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%pri2v',errorFlag)
    END IF ! errorFlag  
   
    WRITE(STDOUT,'(3X,A)') 'Prisms...'    
    READ(iFile) ((gridCENTAUR%pri2v(i,j),j=1,gridCENTAUR%nPris),i=1,6)
    READ(iFile) (idum,i=1,gridCENTAUR%nPris)    
  END IF ! nPrisms

  READ(iFile) gridCENTAUR%nPyrs
  IF ( gridCENTAUR%nPyrs > 0 ) THEN 
    ALLOCATE(gridCENTAUR%pyr2v(5,gridCENTAUR%nPyrs),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%pyr2v',errorFlag)
    END IF ! errorFlag    
  
    WRITE(STDOUT,'(3X,A)') 'Pyramids...'    
    READ(iFile) ((gridCENTAUR%pyr2v(i,j),j=1,gridCENTAUR%nPyrs),i=1,5)
    READ(iFile) (idum,i=1,gridCENTAUR%nPyrs)
  END IF ! nPyramids

  gridCENTAUR%nCells = gridCENTAUR%nTets + gridCENTAUR%nHexs &
                     + gridCENTAUR%nPris + gridCENTAUR%nPyrs

! ==============================================================================
! Boundary types
! ==============================================================================
  
  WRITE(STDOUT,'(3X,A)') 'Boundary information...'   
  
  READ(ifile) gridCENTAUR%nBounds
  
  ALLOCATE(gridCENTAUR%bInfo(3,gridCENTAUR%nBounds),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%bInfo',errorFlag)
  END IF ! errorFlag  

  ALLOCATE(gridCENTAUR%bName(gridCENTAUR%nBounds),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%bName',errorFlag)
  END IF ! errorFlag
    
  READ(ifile) ((gridCENTAUR%bInfo(i,j),j=1,gridCENTAUR%nBounds),i=1,3)

  READ(ifile) (gridCENTAUR%bName(i),i=1,gridCENTAUR%nBounds) 

! ==============================================================================
! Boundary face connectivity
! ==============================================================================

  READ(iFile) gridCENTAUR%nBTris

  IF ( gridCENTAUR%nBTris > 0 ) THEN 
    ALLOCATE(gridCENTAUR%bTri2v(3,gridCENTAUR%nBTris),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%bTri2v',errorFlag)
    END IF ! errorFlag   
  
    WRITE(STDOUT,'(3X,A)') 'Boundary triangles...'    
    READ(iFile) ((gridCENTAUR%bTri2v(i,j),j=1,gridCENTAUR%nBTris),i=1,3)    
  END IF ! gridCENTAUR

  READ(iFile) gridCENTAUR%nBQuads

  IF ( gridCENTAUR%nBQuads > 0 ) THEN
    ALLOCATE(gridCENTAUR%bQuad2v(4,gridCENTAUR%nBQuads),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%bQuad2v',errorFlag)
    END IF ! errorFlag  
   
    WRITE(STDOUT,'(3X,A)') 'Boundary quadrilaterals...'    
    READ(iFile) ((gridCENTAUR%bQuad2v(i,j),j=1,gridCENTAUR%nBQuads),i=1,4)    
  END IF ! gridCENTAUR

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(iFile,IOSTAT=errorFlag)   
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(FILE_CLOSE_ERROR,iFileName,errorFlag)
  END IF  
   
  WRITE(STDOUT,'(1X,A,/)') 'Grid file written successfully.'
  
! ******************************************************************************
! Print grid statistics
! ******************************************************************************

  WRITE(STDOUT,'(/,1X,A)')     'Grid Statistics:'
  WRITE(STDOUT,'(3X,A,2X,I9)') 'Vertices:       ',gridCENTAUR%nVert
  WRITE(STDOUT,'(3X,A,2X,I9)') 'Cells:          ',gridCENTAUR%nCells
  WRITE(STDOUT,'(5X,A,I9)')    'Tetrahedra:     ',gridCENTAUR%nTets
  WRITE(STDOUT,'(5X,A,I9)')    'Hexahedra:      ',gridCENTAUR%nHexs
  WRITE(STDOUT,'(5X,A,I9)')    'Prisms:         ',gridCENTAUR%nPris
  WRITE(STDOUT,'(5X,A,I9)')    'Pyramids:       ',gridCENTAUR%nPyrs
  WRITE(STDOUT,'(3X,A,2X,I9)') 'Boundaries:     ',gridCENTAUR%nBounds
  WRITE(STDOUT,'(3X,A,2X,I9)') 'Boundary faces: ', &
                gridCENTAUR%bInfo(2,gridCENTAUR%nBounds)    &
              + gridCENTAUR%bInfo(3,gridCENTAUR%nBounds)
  WRITE(STDOUT,'(5X,A,I9)')    'Triangles:      ', &
                gridCENTAUR%bInfo(2,gridCENTAUR%nBounds)
  WRITE(STDOUT,'(5X,A,I9)')    'Quadrilaterals: ', &
                gridCENTAUR%bInfo(3,gridCENTAUR%nBounds)

  WRITE(STDOUT,'(1X,A)') 'Boundary statistics:'
  DO i = 1,gridCENTAUR%nBounds
    IF ( i > 1 ) THEN
      triOffs  = gridCENTAUR%bInfo(2,i-1)
      quadOffs = gridCENTAUR%bInfo(3,i-1)           
    ELSE
      triOffs  = 0
      quadOffs = 0      
    END IF ! i    
  
    WRITE(STDOUT,'(3X,I2,1X,A20,1X,I3,2(1X,I6)))') & 
      i,TRIM(gridCENTAUR%bName(i)), &
      gridCENTAUR%bInfo(1,i), & 
      gridCENTAUR%bInfo(2,i) - triOffs, & 
      gridCENTAUR%bInfo(3,i) - quadOffs                      
  END DO ! i

! ******************************************************************************
! End
! ******************************************************************************
   
END SUBROUTINE readGridCENTAUR
