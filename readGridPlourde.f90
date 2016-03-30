! ******************************************************************************
!
! $Id: readGridPlourde.f90,v 1.1 2005/03/05 17:25:53 haselbac Exp $
!
! Filename: readGridPlourde.F90
!
! Purpose: Read 2d grid file in Frederic Plourde format.
!
! Description: None.
!
! Input: None.
! 
! Output: None.
!
! Notes: 
!   1. This routine assumes that the grid is quadrilateral.
!
! Author: Andreas Haselbacher
!
! Copyright: (c) 2002 by the University of Illinois
!
! RCS Revision history:
!
!   $Log: readGridPlourde.f90,v $
!   Revision 1.1  2005/03/05 17:25:53  haselbac
!   Initial revision
!
!   Revision 1.1  2003/03/07 15:26:34  haselbac
!   Initial revision
!
!
! ******************************************************************************

SUBROUTINE readGridPlourde

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

  INTEGER :: i,ib,ie,iFile,iq,it,j,dummyInt,v1l,v2l,v1g,v2g,v3g,v4g,indx, & 
             v1,v2,v3,v4,nc
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: q2bFlag
  CHARACTER :: choice
  CHARACTER*(MAX_STRING_LEN) :: iFileName
  DOUBLE PRECISION :: crossprod,dotprod,x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4, & 
                      xcc,xec,yec,ycc,sx,sy,ct,st,x2d,x3d,x4d,y2d,y3d,y4d, & 
                      x2dd,x3dd,x4dd,y2dd,y3dd,y4dd,theta,fact

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

  WRITE(STDOUT,'(1X,A)') 'Reading grid file in Plourde format...'
  
! ******************************************************************************
! Read grid file
! ******************************************************************************

  READ(iFile,*) grid%nCells
  READ(iFile,*) grid%nVert

  grid%nTris  = 0
  grid%nQuads = grid%nCells
  
  WRITE(STDOUT,'(3X,A,1X,I6)') 'Number of vertices:',grid%nVert
  WRITE(STDOUT,'(3X,A,1X,I6)') 'Number of quadrilaterals:',grid%nQuads

! ==============================================================================
! Read connectivity
! ==============================================================================
 
  ALLOCATE(grid%quad2v(4,grid%nQuads),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%quad2v',errorFlag)
  END IF ! errorFlag

  ALLOCATE(q2bFlag(4,grid%nQuads),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'q2bFlag',errorFlag)
  END IF ! errorFlag

  WRITE(STDOUT,'(3X,A)') 'Connectivity...'

  DO i = 1,grid%nQuads
    READ(iFile,*) grid%quad2v(1,i),grid%quad2v(2,i), & 
                  grid%quad2v(3,i),grid%quad2v(4,i), & 
                  dummyInt,dummyInt,dummyInt,dummyInt, & 
                  q2bFlag(1,i),q2bFlag(2,i),q2bFlag(3,i),q2bFlag(4,i)
  END DO ! i

! ==============================================================================
! Read coordinates 
! ==============================================================================

  ALLOCATE(grid%xy(2,grid%nVert),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%xyz',errorFlag)
  END IF ! errorFlag
  
  WRITE(STDOUT,'(3X,A)') 'Coordinates...'
   
  DO i = 1,grid%nVert
    READ(iFile,*) j,grid%xy(1,i),grid%xy(2,i)
  END DO ! i

! ****************************************************************************** 
! Set remaining quantities
! ******************************************************************************

  WRITE(STDOUT,'(3X,A)') 'Boundaries...'

  grid%nBounds = 4 

  ALLOCATE(grid%bound(grid%nBounds),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound',errorFlag)
  END IF ! errorFlag

  DO i = 1,grid%nQuads
    DO j = 1,4
      IF ( q2bFlag(j,i) == 1 ) THEN ! Outlet
        grid%bound(1)%nEdges = grid%bound(1)%nEdges + 1
      ELSE IF ( q2bFlag(j,i) == 2 ) THEN ! Symmetry axis
        grid%bound(2)%nEdges = grid%bound(2)%nEdges + 1      
      ELSE IF ( q2bFlag(j,i) == 3 ) THEN ! Solid wall
        grid%bound(3)%nEdges = grid%bound(3)%nEdges + 1
      ELSE IF ( q2bFlag(j,i) == 4 ) THEN ! Injection wall
        grid%bound(4)%nEdges = grid%bound(4)%nEdges + 1
      END IF ! q2bFlag 
    END DO ! j
  END DO ! i
  
! ------------------------------------------------------------------------------
! Boundary 1 - Solid wall
! ------------------------------------------------------------------------------ 
 
  grid%bound(1)%bName = 'Solid wall'
  grid%bound(1)%bType = 300 

  ALLOCATE(grid%bound(1)%e2v(2,grid%bound(1)%nEdges),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound%e2v',errorFlag)
  END IF ! errorFlag

  grid%bound(1)%nEdges = 0

  DO i = 1,grid%nQuads
    DO j = 1,4
      IF ( q2bFlag(j,i) == 1 ) THEN ! Outlet
        grid%bound(1)%nEdges = grid%bound(1)%nEdges + 1

        ie = grid%bound(1)%nEdges

        v1l = j
        v2l = j + 1
        
        IF ( v2l > 4 ) THEN 
          v2l = 1
        END IF ! v2l

        grid%bound(1)%e2v(1,ie) = grid%quad2v(v1l,i)
        grid%bound(1)%e2v(2,ie) = grid%quad2v(v2l,i)
        
        v1g = grid%quad2v(1,i)
        v2g = grid%quad2v(2,i)
        v3g = grid%quad2v(3,i)
        v4g = grid%quad2v(4,i)
        
        x1 = grid%xy(1,v1g)
        y1 = grid%xy(2,v1g)             
        x2 = grid%xy(1,v2g)
        y2 = grid%xy(2,v2g)
        x3 = grid%xy(1,v3g)
        y3 = grid%xy(2,v3g)     
        x4 = grid%xy(1,v4g)
        y4 = grid%xy(2,v4g)     

        xcc = 0.25D0*(x1 + x2 + x3 + x4)
        ycc = 0.25D0*(y1 + y2 + y3 + y4)        
               
        x1 = grid%xy(1,grid%quad2v(v1l,i))
        y1 = grid%xy(2,grid%quad2v(v1l,i))
        x2 = grid%xy(1,grid%quad2v(v2l,i))
        y2 = grid%xy(2,grid%quad2v(v2l,i))       
           
        xec = 0.5D0*(x1 + x2)
        yec = 0.5D0*(y1 + y2)
               
        sx = y2 - y1
        sy = x1 - x2
        
        dotprod = sx*(xec-xcc) + sy*(yec-ycc)
        
        IF ( dotprod < 0.0D0 ) THEN 
          grid%bound(1)%e2v(1,ie) = grid%quad2v(v2l,i)
          grid%bound(1)%e2v(2,ie) = grid%quad2v(v1l,i)
        END IF        
      END IF ! q2bFlag 
    END DO ! j
  END DO ! i

! ------------------------------------------------------------------------------
! Boundary 2 - Injection
! ------------------------------------------------------------------------------ 

  grid%bound(2)%bName = 'Injection'
  grid%bound(2)%bType = 300 
  
  ALLOCATE(grid%bound(2)%e2v(2,grid%bound(2)%nEdges),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound%e2v',errorFlag)
  END IF ! errorFlag

  grid%bound(2)%nEdges = 0

  DO i = 1,grid%nQuads
    DO j = 1,4
      IF ( q2bFlag(j,i) == 2 ) THEN ! Symmetry axis
        grid%bound(2)%nEdges = grid%bound(2)%nEdges + 1

        ie = grid%bound(2)%nEdges

        v1l = j
        v2l = j + 1
        
        IF ( v2l > 4 ) THEN 
          v2l = 1
        END IF ! v2l

        grid%bound(2)%e2v(1,ie) = grid%quad2v(v1l,i)
        grid%bound(2)%e2v(2,ie) = grid%quad2v(v2l,i)   
        
        v1g = grid%quad2v(1,i)
        v2g = grid%quad2v(2,i)
        v3g = grid%quad2v(3,i)
        v4g = grid%quad2v(4,i)
        
        x1 = grid%xy(1,v1g)
        y1 = grid%xy(2,v1g)             
        x2 = grid%xy(1,v2g)
        y2 = grid%xy(2,v2g)
        x3 = grid%xy(1,v3g)
        y3 = grid%xy(2,v3g)     
        x4 = grid%xy(1,v4g)
        y4 = grid%xy(2,v4g)     

        xcc = 0.25D0*(x1 + x2 + x3 + x4)
        ycc = 0.25D0*(y1 + y2 + y3 + y4)        
               
        x1 = grid%xy(1,grid%quad2v(v1l,i))
        y1 = grid%xy(2,grid%quad2v(v1l,i))
        x2 = grid%xy(1,grid%quad2v(v2l,i))
        y2 = grid%xy(2,grid%quad2v(v2l,i))       
           
        xec = 0.5D0*(x1 + x2)
        yec = 0.5D0*(y1 + y2)
               
        sx = y2 - y1
        sy = x1 - x2
        
        dotprod = sx*(xec-xcc) + sy*(yec-ycc)
        
        IF ( dotprod < 0.0D0 ) THEN 
          grid%bound(2)%e2v(1,ie) = grid%quad2v(v2l,i)
          grid%bound(2)%e2v(2,ie) = grid%quad2v(v1l,i)   
        END IF          
            
      END IF ! q2bFlag 
    END DO ! j
  END DO ! i

! ------------------------------------------------------------------------------
! Boundary 3 - Symmetry
! ------------------------------------------------------------------------------ 

  grid%bound(3)%bName = 'Symmetry'
  grid%bound(3)%bType = 300 

  ALLOCATE(grid%bound(3)%e2v(2,grid%bound(3)%nEdges),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound%e2v',errorFlag)
  END IF ! errorFlag

  grid%bound(3)%nEdges = 0

  DO i = 1,grid%nQuads
    DO j = 1,4
      IF ( q2bFlag(j,i) == 3 ) THEN ! Solid wall
        grid%bound(3)%nEdges = grid%bound(3)%nEdges + 1

        ie = grid%bound(3)%nEdges

        v1l = j
        v2l = j + 1
        
        IF ( v2l > 4 ) THEN 
          v2l = 1
        END IF ! v2l

        grid%bound(3)%e2v(1,ie) = grid%quad2v(v1l,i)
        grid%bound(3)%e2v(2,ie) = grid%quad2v(v2l,i)
        
        v1g = grid%quad2v(1,i)
        v2g = grid%quad2v(2,i)
        v3g = grid%quad2v(3,i)
        v4g = grid%quad2v(4,i)
        
        x1 = grid%xy(1,v1g)
        y1 = grid%xy(2,v1g)             
        x2 = grid%xy(1,v2g)
        y2 = grid%xy(2,v2g)
        x3 = grid%xy(1,v3g)
        y3 = grid%xy(2,v3g)     
        x4 = grid%xy(1,v4g)
        y4 = grid%xy(2,v4g)     

        xcc = 0.25D0*(x1 + x2 + x3 + x4)
        ycc = 0.25D0*(y1 + y2 + y3 + y4)        
               
        x1 = grid%xy(1,grid%quad2v(v1l,i))
        y1 = grid%xy(2,grid%quad2v(v1l,i))
        x2 = grid%xy(1,grid%quad2v(v2l,i))
        y2 = grid%xy(2,grid%quad2v(v2l,i))       
           
        xec = 0.5D0*(x1 + x2)
        yec = 0.5D0*(y1 + y2)
               
        sx = y2 - y1
        sy = x1 - x2
        
        dotprod = sx*(xec-xcc) + sy*(yec-ycc)
        
        IF ( dotprod < 0.0D0 ) THEN 
          grid%bound(3)%e2v(1,ie) = grid%quad2v(v2l,i)
          grid%bound(3)%e2v(2,ie) = grid%quad2v(v1l,i)
        END IF          
               
      END IF ! q2bFlag 
    END DO ! j
  END DO ! i

! ------------------------------------------------------------------------------
! Boundary 4 - Outlet
! ------------------------------------------------------------------------------ 

  grid%bound(4)%bName = 'Outlet'
  grid%bound(4)%bType = 200 

  ALLOCATE(grid%bound(4)%e2v(2,grid%bound(4)%nEdges),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound%e2v',errorFlag)
  END IF ! errorFlag

  grid%bound(4)%nEdges = 0

  DO i = 1,grid%nQuads
    DO j = 1,4
      IF ( q2bFlag(j,i) == 4 ) THEN ! Solid wall
        grid%bound(4)%nEdges = grid%bound(4)%nEdges + 1

        ie = grid%bound(4)%nEdges

        v1l = j
        v2l = j + 1
        
        IF ( v2l > 4 ) THEN 
          v2l = 1
        END IF ! v2l

        grid%bound(4)%e2v(1,ie) = grid%quad2v(v1l,i)
        grid%bound(4)%e2v(2,ie) = grid%quad2v(v2l,i) 
        
        v1g = grid%quad2v(1,i)
        v2g = grid%quad2v(2,i)
        v3g = grid%quad2v(3,i)
        v4g = grid%quad2v(4,i)
        
        x1 = grid%xy(1,v1g)
        y1 = grid%xy(2,v1g)             
        x2 = grid%xy(1,v2g)
        y2 = grid%xy(2,v2g)
        x3 = grid%xy(1,v3g)
        y3 = grid%xy(2,v3g)     
        x4 = grid%xy(1,v4g)
        y4 = grid%xy(2,v4g)     

        xcc = 0.25D0*(x1 + x2 + x3 + x4)
        ycc = 0.25D0*(y1 + y2 + y3 + y4)        
               
        x1 = grid%xy(1,grid%quad2v(v1l,i))
        y1 = grid%xy(2,grid%quad2v(v1l,i))
        x2 = grid%xy(1,grid%quad2v(v2l,i))
        y2 = grid%xy(2,grid%quad2v(v2l,i))       
           
        xec = 0.5D0*(x1 + x2)
        yec = 0.5D0*(y1 + y2)
               
        sx = y2 - y1
        sy = x1 - x2
        
        dotprod = sx*(xec-xcc) + sy*(yec-ycc)
        
        IF ( dotprod < 0.0D0 ) THEN 
          grid%bound(4)%e2v(1,ie) = grid%quad2v(v2l,i)
          grid%bound(4)%e2v(2,ie) = grid%quad2v(v1l,i) 
        END IF          
              
      END IF ! q2bFlag 
    END DO ! j
  END DO ! i

! BEGIN DEBUG
!  DO ib = 1,grid%nBounds
!    WRITE(STDOUT,*) '*** ',ib,grid%bound(ib)%nEdges
!    DO ie = 1,grid%bound(ib)%nEdges
!      WRITE(STDOUT,*) grid%bound(ib)%e2v(1,ie),grid%bound(ib)%e2v(2,ie)
!    END DO ! ie
!  END DO ! ib
! END DEBUG

! ******************************************************************************
! Check orientation
! ******************************************************************************

  DO i = 1,grid%nQuads
    v1g = grid%quad2v(1,i)
    v2g = grid%quad2v(2,i)
    v3g = grid%quad2v(3,i)
    v4g = grid%quad2v(4,i)
    
    x1 = grid%xy(1,v1g)
    y1 = grid%xy(2,v1g)         
    x2 = grid%xy(1,v2g)
    y2 = grid%xy(2,v2g)
    x3 = grid%xy(1,v3g)
    y3 = grid%xy(2,v3g) 
    x4 = grid%xy(1,v4g)
    y4 = grid%xy(2,v4g) 
    
    crossprod = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
    
    IF ( crossprod < 0.0 ) THEN 
      WRITE(*,*) 'ERROR 1: ',i,crossprod
    END IF ! crossprod  
    
    crossprod = (x3-x1)*(y4-y1) - (x4-x1)*(y3-y1)
    
    IF ( crossprod < 0.0 ) THEN 
      WRITE(*,*) 'ERROR 2: ',i,crossprod
    END IF ! crossprod        
  END DO ! i

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
  
END SUBROUTINE readGridPlourde
