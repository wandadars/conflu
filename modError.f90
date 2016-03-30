! ******************************************************************************
!
! $Id: modError.f90,v 1.3 2005/03/12 03:19:43 haselbac Exp $
!
! Filename: modError.F90
!
! Purpose: Collect definitions of error conditions.
!
! Description:
!
! Input:
! 
! Output:
!
! Notes:
!
! Author: Andreas Haselbacher
!
! Copyright: (c) 2001-2005 by the University of Illinois
!
! RCS Revision history:
!
!   $Log: modError.f90,v $
!   Revision 1.3  2005/03/12 03:19:43  haselbac
!   Removed tab
!
!   Revision 1.2  2005/03/12 03:13:59  haselbac
!   Added error treatment for reading GAMBIT grids
!
!   Revision 1.1  2005/03/05 17:25:53  haselbac
!   Initial revision
!
!   Revision 1.6  2004/12/27 15:31:45  haselbac
!   Added error condition for PLOT3D conversion
!
!   Revision 1.5  2004/10/28 17:02:05  haselbac
!   Added zero width error treatment
!
!   Revision 1.4  2004/10/28 17:00:40  haselbac
!   Added zero width error treatment
!
!   Revision 1.3  2004/07/17 22:11:32  haselbac
!   Added error condition for boundary selection (CENTAUR 3d)
!
!   Revision 1.2  2004/01/13 03:44:23  haselbac
!   Added NBOUNDS_MAX_ERROR
!
!   Revision 1.1  2003/03/07 15:26:34  haselbac
!   Initial revision
!
! ******************************************************************************

MODULE modError

  IMPLICIT NONE
  
  SAVE
  
! ******************************************************************************  
! Basic definitions 
! ******************************************************************************    
  
  INTEGER :: errorFlag
  
  INTEGER, PARAMETER :: ERROR_INCREMENT = 100
  
  INTEGER, PARAMETER :: NO_ERROR = 0
  
! ==============================================================================  
! File handling errors   
! ==============================================================================  
  
  INTEGER, PARAMETER :: FILE_OPEN_ERROR  = NO_ERROR        + ERROR_INCREMENT
  INTEGER, PARAMETER :: FILE_CLOSE_ERROR = FILE_OPEN_ERROR + 1 
  INTEGER, PARAMETER :: FILE_READ_ERROR  = FILE_OPEN_ERROR + 2
  
! ==============================================================================  
! Memory allocation errors   
! ==============================================================================  
  
  INTEGER, PARAMETER :: ALLOCATE_ERROR   = NO_ERROR        + 2*ERROR_INCREMENT
  INTEGER, PARAMETER :: DEALLOCATE_ERROR = ALLOCATE_ERROR  + 1
  
! ==============================================================================  
! Data structure errors   
! ==============================================================================  
  
  INTEGER, PARAMETER :: E2V_ERROR              = NO_ERROR  + 3*ERROR_INCREMENT
  INTEGER, PARAMETER :: NEGATIVE_CVOLUME       = E2V_ERROR + 1 
  INTEGER, PARAMETER :: VOLSUM_ERROR           = E2V_ERROR + 2
  INTEGER, PARAMETER :: BINSEARCH_ERROR        = E2V_ERROR + 3
  INTEGER, PARAMETER :: NORMALSUM_ERROR        = E2V_ERROR + 4
  INTEGER, PARAMETER :: VERTEX_NUMBER_ERROR    = E2V_ERROR + 5
  INTEGER, PARAMETER :: PATCH_NUMBERING_ERROR  = E2V_ERROR + 6
  INTEGER, PARAMETER :: INCORRECT_BFACES_ERROR = E2V_ERROR + 7
  INTEGER, PARAMETER :: BTYPE_CONVERT_ERROR    = E2V_ERROR + 8
  INTEGER, PARAMETER :: INVALID_FACETYPE_ERROR = E2V_ERROR + 9
  INTEGER, PARAMETER :: INVALID_CELLTYPE_ERROR = E2V_ERROR + 10
  INTEGER, PARAMETER :: NBOUNDS_MAX_ERROR      = E2V_ERROR + 11  
  INTEGER, PARAMETER :: NVERT_ERROR            = E2V_ERROR + 12
  INTEGER, PARAMETER :: NTYPE_INVALID_ERROR    = E2V_ERROR + 13
  INTEGER, PARAMETER :: NDP_INVALID_ERROR      = E2V_ERROR + 14
  
! ==============================================================================
! Grid errors
! ==============================================================================

  INTEGER, PARAMETER :: GRID_COBALT_NDIM_ERROR   = NO_ERROR + 4*ERROR_INCREMENT
  INTEGER, PARAMETER :: GRID_COBALT_NZONES_ERROR = GRID_COBALT_NDIM_ERROR + 1
    
! ==============================================================================
! HDF errors
! ==============================================================================  
  
  INTEGER, PARAMETER :: HDF_ERROR        = NO_ERROR + 5*ERROR_INCREMENT
  
! ==============================================================================
! User input (for namelist treatment)
! ==============================================================================  
  
  INTEGER, PARAMETER :: USER_INPUT_ERROR = NO_ERROR + 6*ERROR_INCREMENT
   
! ==============================================================================
! I/O Error
! ==============================================================================  
  
  INTEGER, PARAMETER :: UNEQUAL_NVERTICES_ERROR = NO_ERROR + 7*ERROR_INCREMENT
       
! ==============================================================================
! Miscellaneous
! ==============================================================================  
  
  INTEGER, PARAMETER :: REACHED_DEFAULT  = NO_ERROR + 8*ERROR_INCREMENT
  INTEGER, PARAMETER :: STRING_INVALID   = REACHED_DEFAULT + 1
  INTEGER, PARAMETER :: INFINITE_LOOP    = REACHED_DEFAULT + 2
  
! ==============================================================================
! Conversion errors
! ==============================================================================  
  
  INTEGER, PARAMETER :: BOUNDARY_SELECT_ERROR = NO_ERROR + 9*ERROR_INCREMENT
  INTEGER, PARAMETER :: ZEROWIDTH_ERROR      = BOUNDARY_SELECT_ERROR + 1
  
! ******************************************************************************
! Subroutines
! ******************************************************************************  
  
  CONTAINS
  
! ******************************************************************************
!
! $Id: modError.f90,v 1.3 2005/03/12 03:19:43 haselbac Exp $
!
! Filename: errorHandling.F90
!
! Purpose: Unified handling of error conditions.
!
! Description:
!
! Input:
!   errorCode   Program-internal error designation
!   errorInfo   Additional information, line file or array name
!   errorNumber System error designation, like IOSTAT value
! 
! Output:
!
! Notes: 
!   1. Note optional arguments.
!
! Author: Andreas Haselbacher
!
! Copyright: (c) 2000 by the University of Illinois
!
! RCS Revision history:
!
!   $Log: modError.f90,v $
!   Revision 1.3  2005/03/12 03:19:43  haselbac
!   Removed tab
!
!   Revision 1.2  2005/03/12 03:13:59  haselbac
!   Added error treatment for reading GAMBIT grids
!
!   Revision 1.1  2005/03/05 17:25:53  haselbac
!   Initial revision
!
!   Revision 1.6  2004/12/27 15:31:45  haselbac
!   Added error condition for PLOT3D conversion
!
!   Revision 1.5  2004/10/28 17:02:05  haselbac
!   Added zero width error treatment
!
!   Revision 1.4  2004/10/28 17:00:40  haselbac
!   Added zero width error treatment
!
!   Revision 1.3  2004/07/17 22:11:32  haselbac
!   Added error condition for boundary selection (CENTAUR 3d)
!
!   Revision 1.2  2004/01/13 03:44:23  haselbac
!   Added NBOUNDS_MAX_ERROR
!
!   Revision 1.1  2003/03/07 15:26:34  haselbac
!   Initial revision
!
!
! ******************************************************************************

   SUBROUTINE errorHandling(errorCode,errorInfo,errorNumber)

    USE modGlobals

    IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

    CHARACTER*(MAX_STRING_LEN) :: errorMessage

 ! ==============================================================================
 ! Parameters
 ! ==============================================================================   

    INTEGER, INTENT(IN) :: errorCode
    INTEGER, INTENT(IN), OPTIONAL :: errorNumber 
    CHARACTER(LEN=*), OPTIONAL :: errorInfo

 ! ******************************************************************************
 ! Construct appropriate error message
 ! ******************************************************************************

    SELECT CASE (errorCode)

 ! ==============================================================================
 !   File handling errors
 ! ==============================================================================

      CASE ( FILE_OPEN_ERROR )
        errorMessage = 'Could not open file: '
      CASE ( FILE_CLOSE_ERROR ) 
        errorMessage = 'Could not close file: '
      CASE ( FILE_READ_ERROR ) 
        errorMessage = 'Could not read file: '      

 ! ==============================================================================
 !   Memory allocation errors
 ! ==============================================================================

      CASE ( ALLOCATE_ERROR ) 
        errorMessage = 'Could not allocate memory for: '
      CASE ( DEALLOCATE_ERROR ) 
        errorMessage = 'Could not deallocate memory for: '

 ! ==============================================================================
 !   Data structure errors
 ! ==============================================================================

      CASE ( E2V_ERROR ) 
        errorMessage = 'Inconsistent entries in edgelist: '        
      CASE ( NEGATIVE_CVOLUME ) 
        errorMessage = 'Negative control volumes found.'   
      CASE ( VOLSUM_ERROR )   
        errorMessage = 'Percentage difference of volume sums exceeded limit.'
      CASE ( BINSEARCH_ERROR ) 
        errorMessage = 'Binary search failed.'
      CASE ( NORMALSUM_ERROR ) 
        errorMessage = 'Normal vectors of primal cell do not sum to zero.'
      CASE ( VERTEX_NUMBER_ERROR ) 
        errorMessage = 'Invalid vertex number in array.'
      CASE ( PATCH_NUMBERING_ERROR ) 
        errorMessage = 'Patch numbering inconsistent.'
      CASE ( INCORRECT_BFACES_ERROR ) 
        errorMessage = 'Computed number of boundary faces inconsistent.'      
      CASE ( BTYPE_CONVERT_ERROR ) 
        errorMessage = 'Cannot convert boundary type.'      
      CASE ( INVALID_FACETYPE_ERROR ) 
        errorMessage = 'Face with degree unequal to 3 or 4 detected.'
      CASE ( INVALID_CELLTYPE_ERROR ) 
        errorMessage = 'Invalid cell type detected.'
      CASE ( NBOUNDS_MAX_ERROR ) 
        errorMessage = 'Maximum allowed number of boundaries exceeded.'        
      CASE ( NVERT_ERROR ) 
        errorMessage = 'Numbers of vertices do not match.'
      CASE ( NTYPE_INVALID_ERROR )
        errorMessage = 'Invalid element type.'
      CASE ( NDP_INVALID_ERROR )
        errorMessage = 'Invalid ndp type.'

 ! ==============================================================================
 !   Grid errors
 ! ==============================================================================

      CASE ( GRID_COBALT_NDIM_ERROR ) 
        errorMessage = 'Number of spatial dimensions not equal to 3.'
      CASE ( GRID_COBALT_NZONES_ERROR ) 
        errorMessage = 'Number of zones not equal to 1.'

 ! ==============================================================================
 !   User input
 ! ==============================================================================

      CASE ( USER_INPUT_ERROR ) 
        errorMessage = 'Invalid user input for:'     

! ==============================================================================
!     I/O
! ==============================================================================

      CASE ( UNEQUAL_NVERTICES_ERROR ) 
        errorMessage = 'Number of vertices in flow and grid files different.'  

! ==============================================================================
!   Miscellaneous
! ==============================================================================

      CASE ( REACHED_DEFAULT ) 
        errorMessage = 'Reached default case in routine: '
      CASE ( STRING_INVALID )
        errorMessage = 'Invalid string type.'             
      CASE ( INFINITE_LOOP )
        errorMessage = 'Reached infinite loop.'
        
! ==============================================================================
!   Conversion 
! ==============================================================================         
        
      CASE ( BOUNDARY_SELECT_ERROR ) 
        errorMessage = 'No boundary can be selected for extrusion.'  
      CASE ( ZEROWIDTH_ERROR ) 
        errorMessage = 'Width must be greater than zero.'
        
    END SELECT

 ! ******************************************************************************
 ! Print message and abort
 ! ******************************************************************************

    WRITE(STDOUT,FMT='(/,1X,A)')  '*** ERROR detected in CONFLU *** '

    IF ( PRESENT(errorInfo) ) THEN 
      WRITE(STDOUT,FMT='(1X,A,/)') TRIM(errorMessage)//' '//TRIM(errorInfo)
    ELSE
      WRITE(STDOUT,FMT='(1X,A,/)') TRIM(errorMessage)
    END IF ! PRESENT  

    IF ( PRESENT(errorNumber) ) THEN   
      WRITE(STDOUT,FMT='(1X,A,I6)') 'Internal system error code: ',errorNumber
    END IF ! PRESENT

    WRITE(STDOUT,FMT='(1X,A,/)')  '*** Aborting *** '

    STOP

  END SUBROUTINE errorHandling
  
  
! ******************************************************************************
! End
! ******************************************************************************  
  
END MODULE modError
 
