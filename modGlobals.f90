! ******************************************************************************
!
! $Id: modGlobals.f90,v 1.1 2005/03/05 17:25:53 haselbac Exp $
!
! Filename: modGlobals.F90
!
! Purpose:
!
! Description:
!
! Input:
! 
! Output:
!
! Author: Andreas Haselbacher
!
! Copyright: (c) 2001 by the University of Illinois
!
! RCS Revision history:
!
!   $Log: modGlobals.f90,v $
!   Revision 1.1  2005/03/05 17:25:53  haselbac
!   Initial revision
!
!   Revision 1.3  2005/02/18 20:35:38  haselbac
!   Added is2d flag for 2d PLOT3D grids
!
!   Revision 1.2  2004/01/13 03:45:08  haselbac
!   Added NBOUNDS_MAX, some clean-up
!
!   Revision 1.1  2003/03/07 15:26:34  haselbac
!   Initial revision
!
!
! ******************************************************************************

MODULE modGlobals

  IMPLICIT NONE
  
! ******************************************************************************
! Unit number specifications
! ******************************************************************************

! ==============================================================================
! Read/write I/O
! ==============================================================================

  INTEGER, PARAMETER :: STDIN  = 5, & 
                        STDOUT = 6
 
! ==============================================================================
! File I/O
! ==============================================================================

  INTEGER, PARAMETER :: FILE_UNIT_GRID_INPUT   = 11, & 
                        FILE_UNIT_GRID_OUTPUT  = 21
 
! ******************************************************************************
! Dimensions 
! ******************************************************************************
  
  INTEGER, PARAMETER :: MAX_STRING_LEN = 80

! ******************************************************************************
! Grid stuff 
! ******************************************************************************

  LOGICAL :: is2d

  INTEGER, PARAMETER :: FACE_TYPE_TRI  = 1, & 
                        FACE_TYPE_QUAD = FACE_TYPE_TRI + 1

  INTEGER, PARAMETER :: CELL_TYPE_TETRAHEDRON = 1, & 
                        CELL_TYPE_HEXAHEDRON  = 2, & 
                        CELL_TYPE_PRISM       = 3, & 
                        CELL_TYPE_PYRAMID     = 4

! ******************************************************************************
! Miscellaneous
! ******************************************************************************

  INTEGER, PARAMETER :: DUMMY_ZERO  = 0, & 
                        NBOUNDS_MAX = 10 

  DOUBLE PRECISION, PARAMETER :: RADIAN_PER_DEGREE = 1.74532925E-2

! ******************************************************************************
! End
! ****************************************************************************** 

END MODULE
