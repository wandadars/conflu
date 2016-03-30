! ******************************************************************************
!
! $Id: writeGridCENTAUR.f90,v 1.1 2005/03/05 17:25:53 haselbac Exp $
!
! Filename: writeGridCENTAUR.F90
!
! Purpose: Write grid file from CENTAUR in binary double-precision format.
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
! Copyright: (c) 2001 by the University of Illinois
!
! RCS Revision history:
!
!   $Log: writeGridCENTAUR.f90,v $
!   Revision 1.1  2005/03/05 17:25:53  haselbac
!   Initial revision
!
!   Revision 1.1  2003/03/07 15:26:34  haselbac
!   Initial revision
!
!
! ******************************************************************************

SUBROUTINE writeGridCENTAUR

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

  INTEGER :: i,idum,iFile,j
  CHARACTER*(MAX_STRING_LEN) :: caseName,iFileName

! ******************************************************************************
! Open file and read title
! ******************************************************************************

  WRITE(STDOUT,'(/,1X,A)') 'Writing CENTAUR grid file...'
  WRITE(STDOUT,'(3X,A)') 'Enter case name:'
  READ(STDIN,*) caseName

  iFile     = FILE_UNIT_GRID_INPUT
  iFileName = TRIM(caseName)//'.hyb.bin'

  OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN",IOSTAT=errorFlag)   
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(FILE_OPEN_ERROR,iFileName,errorFlag)
  END IF ! errorFlag

  WRITE(iFile) gridCENTAUR%title
  
! ==============================================================================
! Coordinates
! ==============================================================================

  WRITE(iFile) gridCENTAUR%nVert
  
  WRITE(STDOUT,'(3X,A)') 'Coordinates...'
   
  WRITE(iFile) ((gridCENTAUR%xyz(i,j),j=1,gridCENTAUR%nVert),i=1,3)  

  idum = 0
  WRITE(iFile) (idum,i=1,gridCENTAUR%nVert)

! ==============================================================================
! Cell connectivity
! ==============================================================================

  WRITE(iFile) gridCENTAUR%nTets
  IF ( gridCENTAUR%nTets > 0 ) THEN 
    WRITE(STDOUT,'(3X,A)') 'Tetrahedra...'
    WRITE(iFile) ((gridCENTAUR%tet2v(i,j),j=1,gridCENTAUR%nTets),i=1,4)
    WRITE(iFile) (idum,i=1,gridCENTAUR%nTets)    
  END IF ! nTetrahedra
  
  WRITE(iFile) gridCENTAUR%nHexs
  IF ( gridCENTAUR%nHexs > 0 ) THEN 
    WRITE(STDOUT,'(3X,A)') 'Hexahedra...' 
    WRITE(iFile) ((gridCENTAUR%hex2v(i,j),j=1,gridCENTAUR%nHexs),i=1,8)
    WRITE(iFile) (idum,i=1,gridCENTAUR%nHexs)   
  END IF ! nHexahedra

  WRITE(iFile) gridCENTAUR%nPris
  IF ( gridCENTAUR%nPris > 0 ) THEN 
    WRITE(STDOUT,'(3X,A)') 'Prisms...'    
    WRITE(iFile) ((gridCENTAUR%pri2v(i,j),j=1,gridCENTAUR%nPris),i=1,6)
    WRITE(iFile) (idum,i=1,gridCENTAUR%nPris)    
  END IF ! nPrisms

  WRITE(iFile) gridCENTAUR%nPyrs
  IF ( gridCENTAUR%nPyrs > 0 ) THEN 
    WRITE(STDOUT,'(3X,A)') 'Pyramids...'    
    WRITE(iFile) ((gridCENTAUR%pyr2v(i,j),j=1,gridCENTAUR%nPyrs),i=1,5)
    WRITE(iFile) (idum,i=1,gridCENTAUR%nPyrs)
  END IF ! nPyramids

! ==============================================================================
! Boundary types
! ==============================================================================
  
  WRITE(STDOUT,'(3X,A)') 'Boundary information...'   
  
  WRITE(ifile) gridCENTAUR%nBounds  
  WRITE(ifile) ((gridCENTAUR%bInfo(i,j),j=1,gridCENTAUR%nBounds),i=1,3)

  WRITE(ifile) (gridCENTAUR%bName(i),i=1,gridCENTAUR%nBounds) 

! ==============================================================================
! Boundary face connectivity
! ==============================================================================

  WRITE(iFile) gridCENTAUR%nBTris

  IF ( gridCENTAUR%nBTris > 0 ) THEN 
    WRITE(STDOUT,'(3X,A)') 'Boundary triangles...'    
    WRITE(iFile) ((gridCENTAUR%bTri2v(i,j),j=1,gridCENTAUR%nBTris),i=1,3)    
  END IF ! gridCENTAUR

  WRITE(iFile) gridCENTAUR%nBQuads

  IF ( gridCENTAUR%nBQuads > 0 ) THEN 
    WRITE(STDOUT,'(3X,A)') 'Boundary quadrilaterals...'    
    WRITE(iFile) ((gridCENTAUR%bQuad2v(i,j),j=1,gridCENTAUR%nBQuads),i=1,4)    
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
! End
! ******************************************************************************
   
END SUBROUTINE writeGridCENTAUR
