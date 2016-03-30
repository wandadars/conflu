! ******************************************************************************
!
! $Id: conflu.f90,v 1.8 2016/03/04 16:53:37 neal Exp $
!
! Filename: conflu.F90
!
! Purpose: Main routine of CONFLU.
!
! Description: None.
!
! Input: N/A.
! 
! Output: N/A.
!
! Notes: None.
!
! Author: Andreas Haselbacher
!
! Copyright: (c) 2001 by the University of Illinois
!
! RCS Revision history:
!
!   $Log: conflu.f90,v $
!   Revision 1.8  2016/03/04 16:53:37  neal
!   Added binary grid read/write support option.
!
!   Revision 1.7  2005/03/12 03:19:10  haselbac
!   Fixed missing delimiter
!
!   Revision 1.6  2005/03/12 03:13:25  haselbac
!   Incremented minor version number
!
!   Revision 1.5  2005/03/10 22:06:19  haselbac
!   Added support for GAMBIT 2d neutral files
!
!   Revision 1.4  2005/03/10 01:34:12  haselbac
!   Added ASCII Centaur output
!
!   Revision 1.3  2005/03/09 03:57:29  haselbac
!   Changed minor version number and date
!
!   Revision 1.2  2005/03/07 00:57:14  haselbac
!   Incremented minor version number and date
!
!   Revision 1.1  2005/03/05 17:25:53  haselbac
!   Initial revision
!
!   Revision 1.10  2005/03/02 22:29:53  haselbac
!   Updated version and date
!
!   Revision 1.9  2005/02/18 20:40:09  haselbac
!   Updated version number and date
!
!   Revision 1.8  2005/02/18 20:32:43  haselbac
!   Modified for 2d PLOT3D grids
!
!   Revision 1.7  2004/12/27 15:26:56  haselbac
!   Added PLOT3D option
!
!   Revision 1.5  2004/07/17 22:20:57  haselbac
!   Changed version number and date
!
!   Revision 1.4  2004/07/17 22:10:27  haselbac
!   Added option for CENTAUR grids
!
!   Revision 1.3  2004/02/05 18:19:42  haselbac
!   Added call to distortVertices for SCREAM grids
!
!   Revision 1.2  2004/01/13 03:42:45  haselbac
!   Added TRIANGLE option
!
!   Revision 1.1  2003/03/07 15:26:34  haselbac
!   Initial revision
!
!
! ******************************************************************************

PROGRAM conflu

  USE modGlobals

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
 
  INTEGER :: gridFileFormat

! ******************************************************************************
! Start
! ******************************************************************************
  
  WRITE(STDOUT,'(//,1X,A)') '******************************************'
  WRITE(STDOUT,'(1X,A,/)')  '                 CONFLU                   '
  WRITE(STDOUT,'(1X,A,/)')  '         v2.8.0 of March 11, 2005         '
  WRITE(STDOUT,'(1X,A,/)')  '           Andreas Haselbacher            '
  WRITE(STDOUT,'(1X,A)')    'Center for Simulation of Advanced Rockets '
  WRITE(STDOUT,'(1X,A)')    'University of Illinois at Urbana-Champaign'
  WRITE(STDOUT,'(1X,A,/)')  '******************************************'

! ******************************************************************************
! Get user input
! ******************************************************************************

  WRITE(STDOUT,'(1X,A)')   'Available filters:'
  WRITE(STDOUT,'(1X,A)')   '  Two-dimensional grid formats:'
  WRITE(STDOUT,'(1X,A)')   '    1 - Romuald Carpentier (Ex-INRIA)'
!  WRITE(STDOUT,'(1X,A)')   '    2 - AURORA (Gary Page)' 
  WRITE(STDOUT,'(1X,A)')   '    3 - SCREAM (Andreas Haselbacher)'
  WRITE(STDOUT,'(1X,A)')   '    4 - Frederic Plourde (ENSMA)'
  WRITE(STDOUT,'(1X,A)')   '    5 - Penelope Leyland (EPFL)'
  WRITE(STDOUT,'(1X,A)')   '    6 - TRIANGLE (Jonathan Shewchuk)' 
  WRITE(STDOUT,'(1X,A)')   '    7 - GAMBIT neutral file (FLUENT)' 
  WRITE(STDOUT,'(1X,A)')   '  Three-dimensional grid formats:'
  WRITE(STDOUT,'(1X,A)')   '   10 - COBALT (Matthew Grismer)'
  WRITE(STDOUT,'(1X,A)')   '   11 - CENTAUR (CentaurSoft)'
  WRITE(STDOUT,'(1X,A)')   '   12 - PLOT3D (Pieter Buning, NASA Langley)'
  WRITE(STDOUT,'(1X,A)')   '   13 - PLOT3D (Binary data)'

 
  WRITE(STDOUT,'(/,1X,A)') 'Enter your choice:'
  READ(STDIN,*) gridFileFormat

! ******************************************************************************
! Read grid and convert it to CENTAUR format
! ******************************************************************************

  IF ( gridFileFormat == 1 ) THEN 
    CALL readGridCarpentier
    CALL generic2CENTAUR
!  ELSE IF ( gridFileFormat == 2 ) THEN 
!    CALL readGridAURORA 
!    CALL AURORA2Generic
!    CALL generic2CENTAUR 
  ELSE IF ( gridFileFormat == 3 ) THEN
    CALL readGridSCREAM 
    CALL distortVertices
    CALL generic2CENTAUR
  ELSE IF ( gridFileFormat == 4 ) THEN 
    CALL readGridPlourde 
    CALL generic2CENTAUR 
  ELSE IF ( gridFileFormat == 5 ) THEN 
    CALL readGridLeyland
    CALL generic2CENTAUR
  ELSE IF ( gridFileFormat == 6 ) THEN 
    CALL readGridTRIANGLE
    CALL generic2CENTAUR    
  ELSE IF ( gridFileFormat == 7 ) THEN 
    CALL readGridGAMBIT2d
    CALL generic2CENTAUR    
  ELSE IF ( gridFileFormat == 10 ) THEN
    CALL readGridCOBALT
    CALL COBALT2CENTAUR
  ELSE IF ( gridFileFormat == 11 ) THEN 
    CALL readGridCENTAUR
    CALL CENTAUR2Generic 
    CALL generic2CENTAUR
  ELSE IF ( gridFileFormat == 12 ) THEN 
    CALL readGridPLOT3D
    
    IF ( is2d .EQV. .TRUE. ) THEN 
      CALL PLOT3D2Generic
      CALL generic2CENTAUR
    ELSE  
      CALL PLOT3D2CENTAUR 
    END IF ! is2d
  ELSE IF ( gridFileFormat == 13 ) THEN 
    CALL readGridPLOT3DBinary
 
    IF ( is2d .EQV. .TRUE. ) THEN 
      CALL PLOT3D2Generic
      CALL generic2CENTAUR
    ELSE 
      CALL PLOT3D2CENTAUR 
    END IF ! is2d

  END IF ! gridFileFormat

! ******************************************************************************
! Write grid in CENTAUR format
! ******************************************************************************

  WRITE(STDOUT,'(/,1X,A)') 'Specify output format:'
  WRITE(STDOUT,'(1X,A)') '  1 - CENTAUR ASCII file'
  WRITE(STDOUT,'(1X,A)') '  2 - CENTAUR binary file'  
  WRITE(STDOUT,'(/,1X,A)') 'Enter your choice:'
  READ(STDIN,*) gridFileFormat
  
  IF ( gridFileFormat == 1 ) THEN 
    CALL writeGridCENTAURASCII
  ELSE 
    CALL writeGridCENTAUR
  END IF ! gridFileFormat
  
  WRITE(STDOUT,'(/)')
  
! ******************************************************************************
! End
! ******************************************************************************
  
END PROGRAM conflu
