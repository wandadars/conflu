! ******************************************************************************
!
! $Id: readGridLeyland.f90,v 1.1 2005/03/05 17:25:53 haselbac Exp $
!
! Filename: readGridLeyland.F90
!
! Purpose: Read 2d grid file in Penelope Leylands format.
!
! Description: None.
!
! Input: None.
! 
! Output: None.
!
! Notes: 
!   1. This routine assumes that the grid is triangular.
!   2. This routine is hard-coded for the airintake geometry.
!   3. The third boundary is split into three pieces.
!
! Author: Andreas Haselbacher
!
! Copyright: (c) 2003 by the University of Illinois
!
! RCS Revision history:
!
!   $Log: readGridLeyland.f90,v $
!   Revision 1.1  2005/03/05 17:25:53  haselbac
!   Initial revision
!
!   Revision 1.2  2004/12/27 16:06:01  haselbac
!   Added use of modSortSearch
!
!   Revision 1.1  2003/03/07 15:26:34  haselbac
!   Initial revision
!
!
! ******************************************************************************

SUBROUTINE readGridLeyland

  USE modError
  USE modGlobals
  USE modGrid
  USE modSortSearch

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
! Local variables
! ==============================================================================

  INTEGER :: cntr,c1,dpCntr,dummyInteger,entry,e1,e2,i,ib,ic,ie,iFile,iloc1,iloc2, &
             iq,it,iv,j,lastnnzlink,link,maxv,maxv2,minv,minv2, &
             nBEdges,nBEdges2,nBVert, &
             nBVertIn,nBVertInt,nBVertOut,nBVertSym,nBVertWall,nEdges, &
             nEdgesEst,start,v1,v2,v3,v4
  INTEGER, DIMENSION(:), ALLOCATABLE :: count,strt,vertFlag
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: e2c,e2c_prov,e2v,e2vBak,e2v_prov
  DOUBLE PRECISION :: dp,nm,nx,ny,xa,xm,x1,x2,x3,ya,ym,y1,y2,y3 
  CHARACTER :: choice
  CHARACTER*(MAX_STRING_LEN) :: dummyString,iFileName

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

  WRITE(STDOUT,'(1X,A)') 'Reading grid file in Leyland format...'
  
! ******************************************************************************
! Read grid file
! ******************************************************************************
  
  READ(iFile,'(A)') dummyString
  READ(iFile,'(A)') dummyString
  READ(iFile,'(A)') dummyString  

  READ(iFile,'(A)') dummyString    
  i = SCAN(dummyString,":",.TRUE.)
  READ(dummyString(i+1:LEN_TRIM(dummyString)),*) grid%nVert

  READ(iFile,'(A)') dummyString    
  i = SCAN(dummyString,":",.TRUE.)
  READ(dummyString(i+1:LEN_TRIM(dummyString)),*) grid%nCells
  
  READ(iFile,'(A)') dummyString  
  
  READ(iFile,'(A)') dummyString    
  i = SCAN(dummyString,":",.TRUE.)
  READ(dummyString(i+1:LEN_TRIM(dummyString)),*) nBVert    
  
  grid%nTris  = grid%nCells
  grid%nQuads = 0
  
  WRITE(STDOUT,'(3X,A,1X,I6)') 'Number of vertices:',grid%nVert
  WRITE(STDOUT,'(3X,A,1X,I6)') 'Number of triangles:',grid%nTris

  DO 
    READ(iFile,'(A)') dummyString
    IF ( dummyString == "NODES_2D@" ) THEN 
      EXIT
    END IF ! dummyString
  END DO ! <empty>

  READ(iFile,'(A)') dummyString 

! ==============================================================================
! Read coordinates 
! ==============================================================================

  ALLOCATE(grid%xy(2,grid%nVert),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%xyz',errorFlag)
  END IF ! errorFlag
  
  WRITE(STDOUT,'(3X,A)') 'Coordinates...'
   
  DO i = 1,grid%nVert 
    READ(iFile,'(2(E15.6))') grid%xy(1,i),grid%xy(2,i)
  END DO ! i 

! ==============================================================================
! Read vertex flags
! ==============================================================================

  ALLOCATE(vertFlag(grid%nVert),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'vertFlag',errorFlag)
  END IF ! errorFlag

  READ(iFile,'(A)') dummyString 
  READ(iFile,'(A)') dummyString 

  WRITE(STDOUT,'(3X,A)') 'Vertex flags...'

  DO i = 1,grid%nVert
    READ(iFile,*) vertFlag(i)
  END DO ! i  

! ==============================================================================
! Read levels info (not needed)
! ==============================================================================

  READ(iFile,'(A)') dummyString 
  READ(iFile,'(A)') dummyString 

  DO i = 1,grid%nVert
    READ(iFile,*) dummyInteger
  END DO ! i 

! ==============================================================================
! Read connectivity
! ==============================================================================

  ALLOCATE(grid%tri2v(3,grid%nTris),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%tri2v',errorFlag)
  END IF ! errorFlag

  READ(iFile,'(A)') dummyString 
  READ(iFile,'(A)') dummyString 

  DO i = 1,grid%nTris
    READ(iFile,*) grid%tri2v(1,i),grid%tri2v(2,i),grid%tri2v(3,i)
  END DO ! i 

! ==============================================================================
! Sort vertices
! ==============================================================================

  nBVertIn   = 0
  nBVertInt  = 0  
  nBVertOut  = 0 
  nBVertSym  = 0 
  nBVertWall = 0    

  DO i = 1,grid%nVert
    SELECT CASE ( vertFlag(i) )
      CASE ( -11 ) 
        nBVertOut  = nBVertOut  + 1
        nBVertWall = nBVertWall + 1 
      CASE ( -10 ) 
        nBVertIn   = nBVertIn   + 1
        nBVertWall = nBVertWall + 1         
      CASE ( -8 ) 
        nBVertIn   = nBVertIn   + 1  
      CASE ( -2 ) 
        nBVertWall = nBVertWall + 1               
      CASE ( 0 ) 
        nBVertInt  = nBVertInt  + 1
      CASE ( 1 ) 
        nBVertSym  = nBVertSym  + 1
      CASE ( 2 ) 
        nBVertWall = nBVertWall + 1
      CASE ( 8 ) 
        nBVertIn   = nBVertIn   + 1      
      CASE ( 9 ) 
        nBVertOut  = nBVertOut  + 1    
      CASE DEFAULT
        CALL errorHandling(REACHED_DEFAULT,'readGridLeyland')
    END SELECT ! vertFlag
  END DO ! i

! ==============================================================================
! Define boundary 
! ==============================================================================

  grid%nBounds = 3

  ALLOCATE(grid%bound(grid%nBounds+2),STAT=errorFlag) ! NOTE +2: split below
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound',errorFlag)
  END IF ! errorFlag

  grid%bound(1)%bName = 'Wall'
  grid%bound(1)%bType = 300 

  grid%bound(2)%bName = 'Outlet'
  grid%bound(2)%bType = 200 

  grid%bound(3)%bName = 'InletFront'
  grid%bound(3)%bType = 100 

! ==============================================================================
! Define boundary vertex lists
! ==============================================================================

  DO i = 1,grid%nBounds
    IF ( i == 1 ) THEN ! wall
      ALLOCATE(grid%bound(i)%v(nBVertWall),STAT=errorFlag)
      IF ( errorFlag /= NO_ERROR ) THEN 
        CALL errorHandling(ALLOCATE_ERROR,'grid%bound%v',errorFlag)
      END IF ! errorFlag      
      
      grid%bound(i)%nVert = nBVertWall
      
      cntr = 0
      
      DO j = 1,grid%nVert
        IF ( vertFlag(j) == -11 .OR. & 
             vertFlag(j) == -10 .OR. & 
             vertFlag(j) == - 2 .OR. & 
             vertFlag(j) ==   2 ) THEN 
          cntr = cntr + 1
          grid%bound(i)%v(cntr) = j
        END IF ! vertFlag 
      END DO ! j
      
      IF ( cntr /= grid%bound(i)%nVert ) THEN       
        WRITE(*,*) 'Error - inconsistency!',cntr,grid%bound(i)%nVert
        STOP
      END IF ! cntr
    ELSE IF ( i == 2 ) THEN ! outflow
      ALLOCATE(grid%bound(i)%v(nBVertOut),STAT=errorFlag)
      IF ( errorFlag /= NO_ERROR ) THEN 
        CALL errorHandling(ALLOCATE_ERROR,'grid%bound%v',errorFlag)
      END IF ! errorFlag      
      
      grid%bound(i)%nVert = nBVertOut     
      
      cntr = 0
      
      DO j = 1,grid%nVert
        IF ( vertFlag(j) == -11 .OR. & 
             vertFlag(j) ==   9 ) THEN 
          cntr = cntr + 1             
          grid%bound(i)%v(cntr) = j
        END IF ! vertFlag 
      END DO ! j    
      
      IF ( cntr /= grid%bound(i)%nVert ) THEN       
        WRITE(*,*) 'Error - inconsistency!',cntr,grid%bound(i)%nVert
        STOP
      END IF ! cntr      
    ELSE IF ( i == 3 ) THEN ! inflow
      ALLOCATE(grid%bound(i)%v(nBVertIn),STAT=errorFlag)
      IF ( errorFlag /= NO_ERROR ) THEN 
        CALL errorHandling(ALLOCATE_ERROR,'grid%bound%v',errorFlag)
      END IF ! errorFlag      
      
      grid%bound(i)%nVert = nBVertIn   
            
      cntr = 0
      
      DO j = 1,grid%nVert
        IF ( vertFlag(j) == -10 .OR. & 
             vertFlag(j) ==  -8 .OR. & 
             vertFlag(j) ==   8 ) THEN 
          cntr = cntr + 1             
          grid%bound(i)%v(cntr) = j
        END IF ! vertFlag 
      END DO ! j 
      
      IF ( cntr /= grid%bound(i)%nVert ) THEN       
        WRITE(*,*) 'Error - inconsistency!',cntr,grid%bound(i)%nVert
        STOP
      END IF ! cntr              
    END IF ! i
  END DO ! i

! ******************************************************************************
! Build edge list (from SCREAM)
! ******************************************************************************

  nEdgesEst = 2*3*grid%nVert

  ALLOCATE(e2v_prov(2,nEdgesEst),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'e2v_prov',errorFlag)
  END IF ! errorFlag

  ALLOCATE(count(nEdgesEst),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'count',errorFlag)
  END IF ! errorFlag

  ALLOCATE(strt(grid%nVert),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'strt',errorFlag)
  END IF ! errorFlag

  DO i = 1,nEdgesEst
    e2v_prov(1,i) = 0
    e2v_prov(2,i) = 0
    count(i) = 0  
  ENDDO ! i

  DO i = 1,grid%nVert
    strt(i) = 0
  ENDDO ! i

  nEdges = 0 

! ==============================================================================
! Provisional edge list
! ==============================================================================

  DO ic = 1,grid%nTris
    DO i = 1,3
      e1 = i
      e2 = i + 1

      IF ( e2 .GT. 3 ) THEN 
        e2 = 1
      ENDIF ! e2  

      v1 = grid%tri2v(e1,ic)
      v2 = grid%tri2v(e2,ic)

      IF ( v1 .NE. v2 ) THEN ! to guard against last entry for triangles
        minv = MIN(v1,v2)
        maxv = MAX(v1,v2)

        IF ( strt(minv) .EQ. 0 ) THEN ! first edge at vertex 
          nEdges = nEdges + 1
          strt(minv) = nEdges
          e2v_prov(1,nEdges) = maxv
        ELSE ! other edges exist already 
          lastnnzlink = 0

          entry = e2v_prov(1,strt(minv))
          link  = e2v_prov(2,strt(minv))

          IF ( link .NE. 0 ) THEN 
            lastnnzlink = link
          ENDIF ! link

          DO WHILE ( entry .NE. maxv .AND. link .NE. 0 ) 
            entry = e2v_prov(1,link)
            link  = e2v_prov(2,link)

            IF ( link .NE. 0 ) THEN 
              lastnnzlink = link
            ENDIF ! link
          ENDDO ! while

          IF ( link .EQ. 0 ) THEN 
            IF ( lastnnzlink .EQ. 0 ) THEN 
              link = strt(minv)
            ELSE
              link = lastnnzlink
            ENDIF ! lastnnzlink
          ENDIF ! link

          IF ( entry .NE. maxv ) THEN 
            nEdges = nEdges + 1
            e2v_prov(2,link) = nEdges
            e2v_prov(1,nEdges) = maxv
          ENDIF ! entry
        ENDIF ! strt
      ENDIF ! e1
    ENDDO ! i
  ENDDO ! ic

  ALLOCATE(e2v(2,nEdges),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'e2v',errorFlag)
  END IF ! errorFlag
  
  cntr = 0

  DO iv = 1,grid%nVert
    start = strt(iv)

    IF ( start .NE. 0 ) THEN
      entry = e2v_prov(1,start)
      link  = e2v_prov(2,start)

      DO WHILE ( entry .NE. 0 )
        cntr = cntr + 1

        e2v(1,cntr) = iv
        e2v(2,cntr) = entry

        IF ( link .NE. 0 ) THEN
          entry   = e2v_prov(1,link)
          link    = e2v_prov(2,link)
        ELSE
          entry = 0
        ENDIF ! link
      ENDDO ! while
    ENDIF ! start
  ENDDO ! iv

  DEALLOCATE(e2v_prov,STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(DEALLOCATE_ERROR,'e2v_prov',errorFlag)
  END IF ! errorFlag

  DEALLOCATE(count,STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(DEALLOCATE_ERROR,'count',errorFlag)
  END IF ! errorFlag

  DEALLOCATE(strt,STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(DEALLOCATE_ERROR,'strt',errorFlag)
  END IF ! errorFlag
  
! ==============================================================================
! Build e2c array - very primitive, but works
! ==============================================================================  

  ALLOCATE(e2c(2,nEdges),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'e2c',errorFlag)
  END IF ! errorFlag

  e2c(:,:) = 0

  ALLOCATE(count(nEdges),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'count',errorFlag)
  END IF ! errorFlag

  count(:) = 0

  DO ic = 1,grid%nTris
    DO i = 1,3
      e1 = i
      e2 = i + 1

      IF ( e2 .GT. 3 ) THEN
        e2 = 1
      ENDIF ! e2

      v1 = grid%tri2v(e1,ic)
      v2 = grid%tri2v(e2,ic)

      minv = MIN(v1,v2)
      maxv = MAX(v1,v2)

      DO ie = 1,nEdges
        v3 = e2v(1,ie)
        v4 = e2v(2,ie)

        minv2 = MIN(v3,v4)
        maxv2 = MAX(v3,v4)

        IF ( minv .EQ. minv2 .AND. maxv .EQ. maxv2 ) THEN
          IF ( count(ie) == 2 ) THEN 
            EXIT
          END IF ! count        
        
          count(ie) = count(ie) + 1
          
          IF ( count(ie) <= 2 ) THEN
            e2c(count(ie),ie) = ic
          ELSE
            WRITE(*,*) 'error! count gt 2.'
          ENDIF ! count
          
          IF ( count(ie) == 2 ) THEN 
            EXIT
          END IF ! count          
        ENDIF ! minv

      ENDDO ! ie
    ENDDO ! i
  ENDDO ! ic

  nBEdges = 0

  DO ie = 1,nEdges
    IF ( e2c(2,ie) == 0 ) THEN 
      nBEdges = nBEdges + 1
    END IF ! e2c
  END DO ! ie

  DEALLOCATE(count,STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(DEALLOCATE_ERROR,'count',errorFlag)
  END IF ! errorFlag

! ==============================================================================
! Sort vertices on each boundary, allocate memory for e2v on each boundary
! ==============================================================================

  DO i = 1,grid%nBounds
    CALL quickSortInteger(grid%bound(i)%v,grid%bound(i)%nVert)
  END DO ! i

  DO i = 1,grid%nBounds ! Note upper limit
    ALLOCATE(grid%bound(i)%e2v(2,grid%bound(i)%nVert),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'grid%bound%e2v',errorFlag)
    END IF ! errorFlag
    
    grid%bound(i)%nEdges = 0
    
    ALLOCATE(grid%bound(i)%e2c(grid%bound(i)%nVert),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'grid%bound%e2c',errorFlag)
    END IF ! errorFlag    
  END DO ! i

! ==============================================================================
! Loop over edges and find edges which are on that boundary
! ==============================================================================

  DO ie = 1,nEdges
    v1 = e2v(1,ie)
    v2 = e2v(2,ie)
    
    DO ib = 1,grid%nBounds
      CALL binarySearchInteger(grid%bound(ib)%v,grid%bound(ib)%nVert, &
                               v1,iloc1) 
      
      IF ( iloc1 /= ELEMENT_NOT_FOUND ) THEN 
        CALL binarySearchInteger(grid%bound(ib)%v,grid%bound(ib)%nVert, &
                                 v2,iloc2) 
                                   
        IF ( iloc2 /= ELEMENT_NOT_FOUND ) THEN 
          grid%bound(ib)%nEdges = grid%bound(ib)%nEdges + 1
           
          grid%bound(ib)%e2v(1,grid%bound(ib)%nEdges) = v1
          grid%bound(ib)%e2v(2,grid%bound(ib)%nEdges) = v2
          
          grid%bound(ib)%e2c(grid%bound(ib)%nEdges) = e2c(1,ie)
          
          IF ( e2c(2,ie) /= 0 ) THEN 
            WRITE(*,*) 'Error - e2c(2,ie) should be zero'
            STOP
          END IF ! e2c                
        END IF ! iloc2                                   
      END IF ! iloc1      
    END DO ! ib
  END DO ! ie

  nBEdges2 = 0

  DO i = 1,grid%nBounds
    nBEdges2 = nBEdges2 + grid%bound(i)%nEdges
  END DO ! i 
  
  IF ( nBEdges /= nBEdges2 ) THEN 
    WRITE(*,*) 'Error - inconsistency in number of boundary edges!'
    STOP
  END IF ! nBEdges 

! ******************************************************************************
! Get correct orientation: Normals must point outwards
! ******************************************************************************

  dpCntr = 0 

  DO ib = 1,grid%nBounds
    DO ie = 1,grid%bound(ib)%nEdges
      v1 = grid%bound(ib)%e2v(1,ie)
      v2 = grid%bound(ib)%e2v(2,ie)      
      c1 = grid%bound(ib)%e2c(  ie)
      
      x1 = grid%xy(1,v1)
      y1 = grid%xy(2,v1)
      x2 = grid%xy(1,v2)
      y2 = grid%xy(2,v2)
            
      nx = y2 - y1
      ny = x1 - x2
      
      xm = (x1 + x2)/2.0D0
      ym = (y1 + y2)/2.0D0
      
      v1 = grid%tri2v(1,c1)
      v2 = grid%tri2v(2,c1)
      v3 = grid%tri2v(3,c1)  
      
      x1 = grid%xy(1,v1)
      y1 = grid%xy(2,v1)
      x2 = grid%xy(1,v2)
      y2 = grid%xy(2,v2)
      x3 = grid%xy(1,v3)
      y3 = grid%xy(2,v3)
      
      xa = (x1 + x2 + x3)/3.0D0
      ya = (y1 + y2 + y3)/3.0D0 
      
      dp = (xa-xm)*nx + (ya-ym)*ny
            
      IF ( dp > 0.0D0 ) THEN ! swap
        v1 = grid%bound(ib)%e2v(1,ie)
        v2 = grid%bound(ib)%e2v(2,ie)         
      
        grid%bound(ib)%e2v(1,ie) = v2
        grid%bound(ib)%e2v(2,ie) = v1
        
        dpCntr = dpCntr + 1   
      END IF ! dp                                     
    END DO ! ie
  END DO ! ib

  DO i = 1,grid%nBounds ! Note upper limit
    DEALLOCATE(grid%bound(i)%e2c,STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(DEALLOCATE_ERROR,'grid%bound%e2c',errorFlag)
    END IF ! errorFlag
  END DO ! i

! ******************************************************************************
! Split inflow boundary into three pieces
! ******************************************************************************

  ib = 3
  
  ALLOCATE(e2vBak(2,grid%bound(ib)%nEdges),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'e2vBak',errorFlag)
  END IF ! errorFlag
  
  nEdges = grid%bound(ib)%nEdges
  
  e2vBak(1:2,:nEdges) = grid%bound(ib)%e2v(1:2,1:nEdges)
  
  grid%bound(ib)%e2v(:,:) = 0
    
  grid%bound(ib  )%nEdges = 0
  grid%bound(ib+1)%nEdges = 0
  grid%bound(ib+2)%nEdges = 0
  
  DO ie = 1,nEdges
    v1 = e2vBak(1,ie)
    v2 = e2vBak(2,ie)
  
    x1 = grid%xy(1,v1)
    y1 = grid%xy(2,v1)
    x2 = grid%xy(1,v2)
    y2 = grid%xy(2,v2)

    nx = y2 - y1
    ny = x1 - x2
    nm = SQRT(nx*nx + ny*ny)    
     
    nx = nx/nm
    ny = ny/nm  
      
    IF ( nx < -0.99 ) THEN 
      grid%bound(ib  )%nEdges = grid%bound(ib  )%nEdges + 1
    ELSE IF ( ny > 0.99 ) THEN 
      grid%bound(ib+1)%nEdges = grid%bound(ib+1)%nEdges + 1
    ELSE IF ( nx > 0.99 ) THEN 
      grid%bound(ib+2)%nEdges = grid%bound(ib+2)%nEdges + 1            
    END IF ! 
  END DO ! ie
    
  DEALLOCATE(grid%bound(ib)%e2v,STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(DEALLOCATE_ERROR,'grid%bound(ib)%e2v',errorFlag)
  END IF ! errorFlag
  
  ALLOCATE(grid%bound(ib)%e2v(2,grid%bound(ib)%nEdges),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound(ib)%e2v',errorFlag)
  END IF ! errorFlag    

  ALLOCATE(grid%bound(ib+1)%e2v(2,grid%bound(ib+1)%nEdges),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound(ib)%e2v',errorFlag)
  END IF ! errorFlag 

  ALLOCATE(grid%bound(ib+2)%e2v(2,grid%bound(ib+2)%nEdges),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound(ib)%e2v',errorFlag)
  END IF ! errorFlag 

  grid%bound(ib  )%nEdges = 0
  grid%bound(ib+1)%nEdges = 0
  grid%bound(ib+2)%nEdges = 0

  DO ie = 1,nEdges
    v1 = e2vBak(1,ie)
    v2 = e2vBak(2,ie)
  
    x1 = grid%xy(1,v1)
    y1 = grid%xy(2,v1)
    x2 = grid%xy(1,v2)
    y2 = grid%xy(2,v2)

    nx = y2 - y1
    ny = x1 - x2    
      
    nm = SQRT(nx*nx + ny*ny)    
     
    nx = nx/nm
    ny = ny/nm  
            
    IF ( nx < -0.99 ) THEN 
      grid%bound(ib  )%nEdges = grid%bound(ib  )%nEdges + 1
      grid%bound(ib  )%e2v(1,grid%bound(ib)%nEdges) = v1
      grid%bound(ib  )%e2v(2,grid%bound(ib)%nEdges) = v2      
    ELSE IF ( ny > 0.99 ) THEN 
      grid%bound(ib+1)%nEdges = grid%bound(ib+1)%nEdges + 1
      grid%bound(ib+1)%e2v(1,grid%bound(ib+1)%nEdges) = v1
      grid%bound(ib+1)%e2v(2,grid%bound(ib+1)%nEdges) = v2            
    ELSE IF ( nx > 0.99 ) THEN 
      grid%bound(ib+2)%nEdges = grid%bound(ib+2)%nEdges + 1
      grid%bound(ib+2)%e2v(1,grid%bound(ib+2)%nEdges) = v1
      grid%bound(ib+2)%e2v(2,grid%bound(ib+2)%nEdges) = v2                   
    END IF ! nx
  END DO ! ie

  grid%nBounds = grid%nBounds + 2

  grid%bound(4)%bName = 'InletTop'
  grid%bound(4)%bType = 100 

  grid%bound(5)%bName = 'InletBack'
  grid%bound(5)%bType = 100

! ******************************************************************************
! Clean up
! ******************************************************************************

  DEALLOCATE(e2vBak,STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(DEALLOCATE_ERROR,'e2vBak',errorFlag)
  END IF ! errorFlag

  DEALLOCATE(e2v,STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(DEALLOCATE_ERROR,'e2v',errorFlag)
  END IF ! errorFlag

  DEALLOCATE(e2c,STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(DEALLOCATE_ERROR,'e2c',errorFlag)
  END IF ! errorFlag

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
  
END SUBROUTINE readGridLeyland
