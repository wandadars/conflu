! ******************************************************************************
!
! $Id: writeGridCENTAURASCII.f90,v 1.1 2005/03/10 01:33:39 haselbac Exp $
!
! Filename: writeGridCENTAURASCII.F90
!
! Purpose: Write grid file from CENTAUR in ASCII double-precision format.
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
! Copyright: (c) 2005 by the University of Illinois
!
! RCS Revision history:
!
!   $Log: writeGridCENTAURASCII.f90,v $
!   Revision 1.1  2005/03/10 01:33:39  haselbac
!   Initial revision
!
! ******************************************************************************

SUBROUTINE writeGridCENTAURASCII

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
  iFileName = TRIM(caseName)//'.hyb.asc'

  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN",IOSTAT=errorFlag)   
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(FILE_OPEN_ERROR,iFileName,errorFlag)
  END IF ! errorFlag

  WRITE(iFile,'(A80)') TRIM(gridCENTAUR%title)
  
! ==============================================================================
! Coordinates
! ==============================================================================

  WRITE(iFile,'(I16)') gridCENTAUR%nVert
  
  WRITE(STDOUT,'(3X,A)') 'Coordinates...'
   
  DO i = 1,3   
    WRITE(iFile,'(5E16.9)') (gridCENTAUR%xyz(i,j),j=1,gridCENTAUR%nVert) 
  END DO ! i 

  idum = 0
  WRITE(iFile,'(10I16)') (idum,i=1,gridCENTAUR%nVert)

! ==============================================================================
! Cell connectivity
! ==============================================================================

  WRITE(iFile,'(I16)') gridCENTAUR%nTets
  IF ( gridCENTAUR%nTets > 0 ) THEN 
    WRITE(STDOUT,'(3X,A)') 'Tetrahedra...'
    
    DO i = 1,4
      WRITE(iFile,'(10I16)') (gridCENTAUR%tet2v(i,j),j=1,gridCENTAUR%nTets)
    END DO ! i    
    
    WRITE(iFile,'(10I16)') (idum,i=1,gridCENTAUR%nTets)    
  END IF ! nTetrahedra
  
  WRITE(iFile,'(I16)') gridCENTAUR%nHexs
  IF ( gridCENTAUR%nHexs > 0 ) THEN 
    WRITE(STDOUT,'(3X,A)') 'Hexahedra...' 
    
    DO i = 1,8
      WRITE(iFile,'(10I16)') (gridCENTAUR%hex2v(i,j),j=1,gridCENTAUR%nHexs)
    END DO ! i 
    
    WRITE(iFile,'(10I16)') (idum,i=1,gridCENTAUR%nHexs)   
  END IF ! nHexahedra

  WRITE(iFile,'(I16)') gridCENTAUR%nPris
  IF ( gridCENTAUR%nPris > 0 ) THEN 
    WRITE(STDOUT,'(3X,A)') 'Prisms...'    
    
    DO i = 1,6
      WRITE(iFile,'(10I16)') (gridCENTAUR%pri2v(i,j),j=1,gridCENTAUR%nPris)
    END DO ! i
 
    WRITE(iFile,'(10I16)') (idum,i=1,gridCENTAUR%nPris)    
  END IF ! nPrisms

  WRITE(iFile,'(I16)') gridCENTAUR%nPyrs
  IF ( gridCENTAUR%nPyrs > 0 ) THEN 
    WRITE(STDOUT,'(3X,A)') 'Pyramids...'  
    
    DO i = 1,5  
      WRITE(iFile,'(10I16)') (gridCENTAUR%pyr2v(i,j),j=1,gridCENTAUR%nPyrs)
    END DO ! i 
    
    WRITE(iFile,'(10I16)') (idum,i=1,gridCENTAUR%nPyrs)
  END IF ! nPyramids

! ==============================================================================
! Boundary types
! ==============================================================================
  
  WRITE(STDOUT,'(3X,A)') 'Boundary information...'   
  
  WRITE(ifile,'(I16)') gridCENTAUR%nBounds  
  
  DO i = 1,3
    WRITE(ifile,'(10I16)') (gridCENTAUR%bInfo(i,j),j=1,gridCENTAUR%nBounds)
  END DO ! i

  WRITE(ifile,'(A80)') (gridCENTAUR%bName(i),i=1,gridCENTAUR%nBounds) 

! ==============================================================================
! Boundary face connectivity
! ==============================================================================

  WRITE(iFile,'(I16)') gridCENTAUR%nBTris

  IF ( gridCENTAUR%nBTris > 0 ) THEN 
    WRITE(STDOUT,'(3X,A)') 'Boundary triangles...'    
    
    DO i = 1,3
      WRITE(iFile,'(10I16)') (gridCENTAUR%bTri2v(i,j),j=1,gridCENTAUR%nBTris)    
    END DO ! i
  END IF ! gridCENTAUR

  WRITE(iFile,'(I16)') gridCENTAUR%nBQuads

  IF ( gridCENTAUR%nBQuads > 0 ) THEN 
    WRITE(STDOUT,'(3X,A)') 'Boundary quadrilaterals...'    
    
    DO i = 1,4
      WRITE(iFile,'(10I16)') (gridCENTAUR%bQuad2v(i,j),j=1,gridCENTAUR%nBQuads)
    END DO ! i
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
   
END SUBROUTINE writeGridCENTAURASCII
