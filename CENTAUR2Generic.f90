! ******************************************************************************
!
! $Id: CENTAUR2Generic.f90,v 1.2 2005/03/07 00:55:09 haselbac Exp $
!
! Filename: CENTAUR2Generic.F90
!
! Purpose: Convert 3d grid in CENTAUR format into 2d grid.
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
!   $Log: CENTAUR2Generic.f90,v $
!   Revision 1.2  2005/03/07 00:55:09  haselbac
!   Added legend to planarity output
!
!   Revision 1.1  2005/03/05 17:25:53  haselbac
!   Initial revision
!
!   Revision 1.3  2004/12/27 16:05:44  haselbac
!   Added use of modSortSearch
!
!   Revision 1.2  2004/07/19 18:54:44  haselbac
!   Added deletion of empty patches and mapping of patches
!
!   Revision 1.1  2004/07/17 22:17:40  haselbac
!   Initial revision
!
!
! ******************************************************************************

SUBROUTINE CENTAUR2Generic

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
 
  LOGICAL :: xFlatFlag,yFlatFlag,zFlatFlag
  INTEGER :: cntr,c1,c2,entry,entryLoc,e1,e2,flipFlag,i,ib,ib2,ib2d,ib3d, & 
             ic,ie,ie2,iq,iqBeg,iqEnd,it,itBeg,itEnd,iv,ivg,j,lastnnzlink, & 
             link,mapFlag,maxv,minv,quadOffs,nBounds,nBQuads,nBTris,nEdges, & 
             nEdgesEst,nVertEst,selectFlagCntr, &
             start,surfNormDir,term,triOffs,v,v1,v1g,v12,v12g,v2,v2g,v22, &
             v22g,v3,v3g,v4,v4g,v5,v5g,v6,v6g
  INTEGER, DIMENSION(:), ALLOCATABLE :: bCntr,boundMap,flag,flip,strt
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: bPntr,e2c,e2v
  DOUBLE PRECISION :: nmag,nx,nxmax,nxmin,ny,nymax,nymin,nz,nzmax,nzmin, & 
                      toler,x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,z1,z2,z3,z4

! ******************************************************************************
! Start, get additional information for extrusion
! ******************************************************************************

  WRITE(STDOUT,'(/,1X,A)') 'Creating 2d-grid in generic format...'

! ******************************************************************************
! Check for empty patches
! ******************************************************************************

  WRITE(STDOUT,'(3X,A)') 'Checking for empty patches...'

  grid3d%nBounds = gridCENTAUR%nBounds
  nBounds        = gridCENTAUR%nBounds

  ib = 0

  IF ( grid3d%nBounds > 0 ) THEN 
    ib = ib + 1

    DO ib2 = 1,nBounds
      IF ( ib2 /= 1 ) THEN 
        triOffs  = gridCENTAUR%bInfo(2,ib2-1) 
        quadOffs = gridCENTAUR%bInfo(3,ib2-1) 
      ELSE 
        triOffs  = 0
        quadOffs = 0
      END IF ! ib2

      nBTris  = gridCENTAUR%bInfo(2,ib2) - triOffs
      nBQuads = gridCENTAUR%bInfo(3,ib2) - quadOffs             

      IF ( (nBTris + nBQuads) /= 0 ) THEN
        gridCENTAUR%bInfo(1,ib) = gridCENTAUR%bInfo(1,ib2)
        gridCENTAUR%bInfo(2,ib) = gridCENTAUR%bInfo(2,ib2)
        gridCENTAUR%bInfo(3,ib) = gridCENTAUR%bInfo(3,ib2)

        ib = ib + 1      
      ELSE 
        WRITE(STDOUT,'(3X,A,1X,I3,1X,A)') & 
          '*** WARNING *** Patch',ib2,'is empty and will be deleted!' 

        grid3d%nBounds = grid3d%nBounds - 1     
      END IF ! nBTris    
    END DO ! ib2

    gridCENTAUR%nBounds = grid3d%nBounds
  END IF ! grid3d%nBounds  

! ******************************************************************************
! Build mapping from CENTAUR patches to new patches
! ******************************************************************************  

  ALLOCATE(boundMap(grid3d%nBounds),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'boundMap',errorFlag)
  END IF ! errorFlag

  WRITE(STDOUT,'(3X,A)') 'Would you like to map patches? (0/1)'
  READ(STDIN,*) mapFlag

  IF ( mapFlag == 0 ) THEN 
    DO ib = 1,grid3d%nBounds 
      boundMap(ib) = ib
    END DO ! ib
  ELSE 
    ib = 1
  
    emptyLoop: DO 
      WRITE(STDOUT,'(5X,A,1X,I2,1X,A)') 'Enter new patch number for patch:',ib
      READ(STDIN,*) boundMap(ib)
      
      IF ( boundMap(ib) < 0 .OR. boundMap(ib) > grid3d%nBounds ) THEN 
        WRITE(STDOUT,*) 'Invalid input. Try again.'
      ELSE 
        ib = ib + 1        
      END IF ! boundMap
      
      IF ( ib > grid3d%nBounds ) THEN 
        EXIT emptyLoop
      END IF ! ib
    END DO emptyLoop
    
    ALLOCATE(bCntr(grid3d%nBounds),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'bCntr',errorFlag)
    END IF ! errorFlag  
    
    DO ib = 1,grid3d%nBounds
      bCntr(ib) = 0         
    END DO ! ib      
    
    DO ib = 1,grid3d%nBounds
      bCntr(boundMap(ib)) = bCntr(boundMap(ib)) + 1 
    END DO ! ib  
        
    nBounds = 0
    
    DO ib = 1,grid3d%nBounds
      IF ( bCntr(ib) /= 0 ) THEN 
        nBounds = nBounds + 1
      END IF ! bCntr         
    END DO ! ib    
    
    grid3d%nBounds = nBounds
    
    WRITE(STDOUT,'(5X,A,1X,I2,1X,A)') 'There are',grid3d%nBounds, & 
                                      'patches after mapping.'
                                      
    DEALLOCATE(bCntr,STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(DEALLOCATE_ERROR,'bCntr',errorFlag)
    END IF ! errorFlag                                                
  END IF ! mapFlag 
  
  WRITE(STDOUT,'(5X,A)') 'Mapping:'
  
  DO ib = 1,gridCENTAUR%nBounds 
    WRITE(STDOUT,'(7X,2(1X,I2))') ib,boundMap(ib)
  END DO ! ib  

! ******************************************************************************
! Allocate memory for 3d grid data structure
! ******************************************************************************  
  
  ALLOCATE(grid3d%bound(grid3d%nBounds),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid3d%bound',errorFlag)
  END IF ! errorFlag  
  
! ******************************************************************************
! Extract connectivity to 3d boundary grid
! ******************************************************************************  
  
  WRITE(STDOUT,'(3X,A)') 'Extracting connectivity...'
  
  DO ib = 1,grid3d%nBounds
    grid3d%bound(ib)%nTris  = 0 
    grid3d%bound(ib)%nQuads = 0        
  END DO ! ib  

! ==============================================================================
! Determine number of faces on each boundary 
! ==============================================================================
  
  DO ib = 1,gridCENTAUR%nBounds
    IF ( ib > 1 ) THEN
      triOffs  = gridCENTAUR%bInfo(2,ib-1)
      quadOffs = gridCENTAUR%bInfo(3,ib-1)           
    ELSE
      triOffs  = 0
      quadOffs = 0      
    END IF ! ib    

    nBTris  = gridCENTAUR%bInfo(2,ib) - triOffs 
    nBQuads = gridCENTAUR%bInfo(3,ib) - quadOffs
    
    ib2 = boundMap(ib)
    
    grid3d%bound(ib2)%nTris  = grid3d%bound(ib2)%nTris  + nBTris
    grid3d%bound(ib2)%nQuads = grid3d%bound(ib2)%nQuads + nBQuads        
  END DO ! ib

  WRITE(STDOUT,'(3X,A)') 'Patch sizes after mapping:'
  
  DO ib = 1,grid3d%nBounds
    WRITE(STDOUT,'(5X,I2,2(1X,I6))') ib,grid3d%bound(ib)%nTris, & 
                                     grid3d%bound(ib)%nQuads        
  END DO ! ib   
  
! ==============================================================================
! Allocate memory 
! ==============================================================================
  
  DO ib = 1,grid3d%nBounds  
    IF ( grid3d%bound(ib)%nTris > 0 ) THEN 
      ALLOCATE(grid3d%bound(ib)%tri2v(3,grid3d%bound(ib)%nTris),STAT=errorFlag)
      IF ( errorFlag /= NO_ERROR ) THEN 
        CALL errorHandling(ALLOCATE_ERROR,'grid3d%bound(ib)%tri2v',errorFlag)
      END IF ! errorFlag
    END IF ! grid3d%bound(ib)%nTris

    IF ( grid3d%bound(ib)%nQuads > 0 ) THEN 
      ALLOCATE(grid3d%bound(ib)%Quad2v(4,grid3d%bound(ib)%nQuads), & 
               STAT=errorFlag)
      IF ( errorFlag /= NO_ERROR ) THEN 
        CALL errorHandling(ALLOCATE_ERROR,'grid3d%bound(ib)%quad2v',errorFlag)
      END IF ! errorFlag
    END IF ! grid3d%bound(ib)%nQuads  
  END DO ! ib
    
! ==============================================================================
! Set connectivity 
! ==============================================================================
       
  DO ib = 1,grid3d%nBounds
    grid3d%bound(ib)%nTris  = 0 
    grid3d%bound(ib)%nQuads = 0        
  END DO ! ib       
        
  DO ib = 1,gridCENTAUR%nBounds
    IF ( ib > 1 ) THEN
      itBeg = gridCENTAUR%bInfo(2,ib-1) + 1
      iqBeg = gridCENTAUR%bInfo(3,ib-1) + 1           
    ELSE
      itBeg = 1
      iqBeg = 1      
    END IF ! ib    

    itEnd = gridCENTAUR%bInfo(2,ib) 
    iqEnd = gridCENTAUR%bInfo(3,ib) 
        
    ib2 = boundMap(ib)

    triOffs  = grid3d%bound(ib2)%nTris  
    quadOffs = grid3d%bound(ib2)%nQuads        
    
    DO it = 1,itEnd-itBeg+1
      grid3d%bound(ib2)%tri2v(1,triOffs+it) = gridCENTAUR%bTri2v(1,itBeg+it-1)
      grid3d%bound(ib2)%tri2v(2,triOffs+it) = gridCENTAUR%bTri2v(2,itBeg+it-1)
      grid3d%bound(ib2)%tri2v(3,triOffs+it) = gridCENTAUR%bTri2v(3,itBeg+it-1)            
    END DO ! it

    DO iq = 1,iqEnd-iqBeg+1
      grid3d%bound(ib2)%quad2v(1,quadOffs+iq) = gridCENTAUR%bQuad2v(1,iqBeg+iq-1)
      grid3d%bound(ib2)%quad2v(2,quadOffs+iq) = gridCENTAUR%bQuad2v(2,iqBeg+iq-1)
      grid3d%bound(ib2)%quad2v(3,quadOffs+iq) = gridCENTAUR%bQuad2v(3,iqBeg+iq-1)
      grid3d%bound(ib2)%quad2v(4,quadOffs+iq) = gridCENTAUR%bQuad2v(4,iqBeg+iq-1)                  
    END DO ! iq
        
    grid3d%bound(ib2)%nTris  = grid3d%bound(ib2)%nTris  + itEnd-itBeg+1
    grid3d%bound(ib2)%nQuads = grid3d%bound(ib2)%nQuads + iqEnd-iqBeg+1 
  END DO ! ib    
      
! ******************************************************************************
! Build vertex list on each patch
! ******************************************************************************  

  WRITE(STDOUT,'(3X,A)') 'Building vertex lists...'
      
  DO ib = 1,grid3d%nBounds 
  
! ==============================================================================
!   Determine size
! ==============================================================================  
   
    nVertEst = gridCENTAUR%nVert

    ALLOCATE(flag(nVertEst),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'flag',errorFlag)
    END IF ! errorFlag
        
    grid3d%bound(ib)%nVert = 0    

    DO iv = 1,nVertEst
      flag(iv) = 0
    END DO ! iv
            
    DO it = 1,grid3d%bound(ib)%nTris
      DO i = 1,3
        v = grid3d%bound(ib)%tri2v(i,it)
        
        IF ( flag(v) == 0 ) THEN 
          flag(v) = 1
          grid3d%bound(ib)%nVert = grid3d%bound(ib)%nVert + 1
        END IF ! vFlag  
      END DO ! i
    END DO ! it    
    
    DO iq = 1,grid3d%bound(ib)%nQuads
      DO i = 1,4
        v = grid3d%bound(ib)%quad2v(i,iq)
        
        IF ( flag(v) == 0 ) THEN 
          flag(v) = 1
          grid3d%bound(ib)%nVert = grid3d%bound(ib)%nVert + 1
        END IF ! vFlag  
      END DO ! i
    END DO ! iq
    
    WRITE(STDOUT,'(5X,I2,1X,I5)') ib,grid3d%bound(ib)%nVert    

! ==============================================================================
!   Determine actual list
! ==============================================================================  

    ALLOCATE(grid3d%bound(ib)%v(grid3d%bound(ib)%nVert),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'grid3d%bound(ib)%v',errorFlag)
    END IF ! errorFlag

    DO iv = 1,nVertEst
      flag(iv) = 0
    END DO ! iv

    grid3d%bound(ib)%nVert = 0

    DO it = 1,grid3d%bound(ib)%nTris
      DO i = 1,3
        v = grid3d%bound(ib)%tri2v(i,it)
        
        IF ( flag(v) == 0 ) THEN 
          flag(v) = 1
          grid3d%bound(ib)%nVert = grid3d%bound(ib)%nVert + 1
          grid3d%bound(ib)%v(grid3d%bound(ib)%nVert) = v
        END IF ! vFlag  
      END DO ! i
    END DO ! it    

    DO iq = 1,grid3d%bound(ib)%nQuads
      DO i = 1,4
        v = grid3d%bound(ib)%quad2v(i,iq)
        
        IF ( flag(v) == 0 ) THEN 
          flag(v) = 1
          grid3d%bound(ib)%nVert = grid3d%bound(ib)%nVert + 1
          grid3d%bound(ib)%v(grid3d%bound(ib)%nVert) = v          
        END IF ! vFlag  
      END DO ! i
    END DO ! iq

    DEALLOCATE(flag,STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(DEALLOCATE_ERROR,'flag',errorFlag)
    END IF ! errorFlag
    
! ==============================================================================
!   Sort into ascending order
! ==============================================================================  
    
    CALL quickSortInteger(grid3d%bound(ib)%v,grid3d%bound(ib)%nVert)
  END DO ! ib      
      
! ******************************************************************************
! Renumber connectivity lists
! ******************************************************************************  

  WRITE(STDOUT,'(3X,A)') 'Renumbering connectivity...'
      
  DO ib = 1,grid3d%nBounds  
    DO it = 1,grid3d%bound(ib)%nTris
      DO i = 1,3
        v = grid3d%bound(ib)%tri2v(i,it)
        
        CALL binarySearchInteger(grid3d%bound(ib)%v,grid3d%bound(ib)%nVert,v,j)
        
        IF ( j /= ELEMENT_NOT_FOUND ) THEN 
          grid3d%bound(ib)%tri2v(i,it) = j
        ELSE  
          CALL errorHandling(BINSEARCH_ERROR)
        END IF ! i  
      END DO ! i
    END DO ! it    

    DO iq = 1,grid3d%bound(ib)%nQuads
      DO i = 1,4
        v = grid3d%bound(ib)%quad2v(i,iq)
                
        CALL binarySearchInteger(grid3d%bound(ib)%v,grid3d%bound(ib)%nVert,v,j)
        
        IF ( j /= ELEMENT_NOT_FOUND ) THEN 
          grid3d%bound(ib)%quad2v(i,iq) = j
        ELSE  
          CALL errorHandling(BINSEARCH_ERROR)
        END IF ! i 
      END DO ! i
    END DO ! iq    
  END DO ! ib      
      
! ******************************************************************************
! Build edge list on each patch. NOTE need to track whether edge vertices are
! flipped (flipFlag) when inserting into linked list because need to make sure
! that retain orientation of edges on perimeter of boundaries. 
! ******************************************************************************  
  
  WRITE(STDOUT,'(3X,A)') 'Building edge lists...'
  
  DO ib = 1,grid3d%nBounds  
    nEdgesEst = grid3d%bound(ib)%nVert & 
              + grid3d%bound(ib)%nTris & 
              + grid3d%bound(ib)%nQuads & 
              + 100 

    ALLOCATE(e2c(2,nEdgesEst),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'e2c',errorFlag)
    END IF ! errorFlag

    DO ie = 1,nEdgesEst
      e2c(1,ie) = 0
      e2c(2,ie) = 0 
    END DO ! ie
    
    ALLOCATE(e2v(2,nEdgesEst),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'e2v',errorFlag)
    END IF ! errorFlag
    
    DO ie = 1,nEdgesEst
      e2v(1,ie) = 0
      e2v(2,ie) = 0 
    END DO ! ie
    
    ALLOCATE(strt(grid3d%bound(ib)%nVert),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'strt',errorFlag)
    END IF ! errorFlag

    ALLOCATE(flip(nEdgesEst),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'flip',errorFlag)
    END IF ! errorFlag
    
    DO iv = 1,grid3d%bound(ib)%nVert
      strt(iv) = 0 
    END DO ! iv
   
    nEdges = 0

! ==============================================================================
!   Triangles
! ==============================================================================  
    
    DO it = 1,grid3d%bound(ib)%nTris
      DO i = 1,3
        e1 = i
        e2 = i + 1

        IF ( e2 .GT. 3 ) THEN 
          e2 = 1
        ENDIF ! e2  

        v1 = grid3d%bound(ib)%tri2v(e1,it)
        v2 = grid3d%bound(ib)%tri2v(e2,it)

        minv = MIN(v1,v2)
        maxv = MAX(v1,v2)

        IF ( v1 /= minv ) THEN 
          flipFlag = 1
        ELSE 
          flipFlag = 0 
        END IF ! v1

        IF ( strt(minv) .EQ. 0 ) THEN ! first edge at vertex 
          nEdges = nEdges + 1                              
          strt(minv) = nEdges
          e2v(1,nEdges) = maxv
          e2c(1,nEdges) = it
          flip(nEdges) = flipFlag                              
        ELSE ! other edges exist already 
          lastnnzlink = 0
          entry = e2v(1,strt(minv))
          link  = e2v(2,strt(minv))
          entryLoc = strt(minv)

          IF ( link .NE. 0 ) THEN 
            lastnnzlink = link
          ENDIF ! link

          DO WHILE ( entry .NE. maxv .AND. link .NE. 0 ) 
            entry    = e2v(1,link)
            entryLoc = link            
            link     = e2v(2,link)

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
            e2v(2,link) = nEdges
            e2v(1,nEdges) = maxv
            e2c(1,nEdges) = it 
            flip(nEdges) = flipFlag                      
          ELSE 
            e2c(2,entryLoc) = it
          ENDIF ! entry 
        ENDIF ! strt
      ENDDO ! i
    ENDDO ! it

! ==============================================================================
!   Quadrilaterals
! ==============================================================================  

    quadOffs = grid3d%bound(ib)%nTris

    DO iq = 1,grid3d%bound(ib)%nQuads
      DO i = 1,4
        e1 = i
        e2 = i + 1

        IF ( e2 .GT. 4 ) THEN 
          e2 = 1
        ENDIF ! e2  

        v1 = grid3d%bound(ib)%quad2v(e1,iq)
        v2 = grid3d%bound(ib)%quad2v(e2,iq)

        minv = MIN(v1,v2)
        maxv = MAX(v1,v2)

        IF ( v1 /= minv ) THEN 
          flipFlag = 1
        ELSE 
          flipFlag = 0 
        END IF ! v1

        IF ( strt(minv) .EQ. 0 ) THEN ! first edge at vertex 
          nEdges = nEdges + 1                              
          strt(minv) = nEdges
          e2v(1,nEdges) = maxv
          e2c(1,nEdges) = iq + quadOffs
          flip(nEdges) = flipFlag            
        ELSE ! other edges exist already 
          lastnnzlink = 0

          entry = e2v(1,strt(minv))
          link  = e2v(2,strt(minv))
          entryLoc = strt(minv)

          IF ( link .NE. 0 ) THEN 
            lastnnzlink = link
          ENDIF ! link

          DO WHILE ( entry .NE. maxv .AND. link .NE. 0 ) 
            entry    = e2v(1,link)
            entryLoc = link                  
            link     = e2v(2,link)

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
            e2v(2,link) = nEdges
            e2v(1,nEdges) = maxv
            e2c(1,nEdges) = iq + quadOffs
            flip(nEdges) = flipFlag                           
          ELSE 
            e2c(2,entryLoc) = iq + quadOffs
          ENDIF ! entry                         
        ENDIF ! strt
      ENDDO ! i
    ENDDO ! iq

    WRITE(STDOUT,'(5X,I2,1X,I5)') ib,nEdges

    grid3d%bound(ib)%nEdges = nEdges

    ALLOCATE(grid3d%bound(ib)%e2v(2,grid3d%bound(ib)%nEdges),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'grid3d%bound(ib)%e2v',errorFlag)
    END IF ! errorFlag    

    ALLOCATE(grid3d%bound(ib)%e2c(2,grid3d%bound(ib)%nEdges),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'grid3d%bound(ib)%e2c',errorFlag)
    END IF ! errorFlag    

    cntr = 0
  
    DO iv = 1,grid3d%bound(ib)%nVert    
      start = strt(iv)

      IF ( start .NE. 0 ) THEN
        entry = e2v(1,start)
        link  = e2v(2,start)

        c1 = e2c(1,start)
        c2 = e2c(2,start)
        
        flipFlag = flip(start)

        DO WHILE ( entry .NE. 0 )
          cntr = cntr + 1

          IF ( flipFlag == 0 ) THEN 
            grid3d%bound(ib)%e2v(1,cntr) = iv
            grid3d%bound(ib)%e2v(2,cntr) = entry
          ELSE 
            grid3d%bound(ib)%e2v(2,cntr) = iv
            grid3d%bound(ib)%e2v(1,cntr) = entry            
          END IF ! flipFlag

          grid3d%bound(ib)%e2c(1,cntr) = c1
          grid3d%bound(ib)%e2c(2,cntr) = c2

          IF ( link .NE. 0 ) THEN
            entry = e2v(1,link)        
            c1 = e2c(1,link)
            c2 = e2c(2,link) 
            flipFlag = flip(link)           
            
            link  = e2v(2,link)
          ELSE
            entry = 0
          ENDIF ! link
        ENDDO ! while
      ENDIF ! start
    ENDDO ! iv
        
    DEALLOCATE(e2c,STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(DEALLOCATE_ERROR,'e2c',errorFlag)
    END IF ! errorFlag

    DEALLOCATE(e2v,STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(DEALLOCATE_ERROR,'e2v',errorFlag)
    END IF ! errorFlag
    
    DEALLOCATE(strt,STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(DEALLOCATE_ERROR,'strt',errorFlag)
    END IF ! errorFlag

    DEALLOCATE(flip,STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(DEALLOCATE_ERROR,'flip',errorFlag)
    END IF ! errorFlag
  END DO ! ib
     
! ******************************************************************************
! Determine number of boundary edges. Build list of edge-to-vertex and edge-to-
! cell connectivity arrays for boundary edges. 
! ******************************************************************************     
   
  WRITE(STDOUT,'(3X,A)') 'Building bounding edge list...'
   
  DO ib = 1,grid3d%nBounds  
    grid3d%bound(ib)%nBEdges = 0
  
    DO ie = 1,grid3d%bound(ib)%nEdges        
      IF ( grid3d%bound(ib)%e2c(1,ie) == 0 .OR. & 
           grid3d%bound(ib)%e2c(2,ie) == 0 ) THEN 
        grid3d%bound(ib)%nBEdges = grid3d%bound(ib)%nBEdges + 1
      ELSE IF ( grid3d%bound(ib)%e2c(1,ie) == 0 .AND. & 
                grid3d%bound(ib)%e2c(2,ie) == 0 ) THEN
        STOP
      END IF ! c1
    END DO ! ie
    
    WRITE(STDOUT,'(5X,I2,1X,I5)') ib,grid3d%bound(ib)%nBEdges

    ALLOCATE(grid3d%bound(ib)%be2v(2,grid3d%bound(ib)%nBEdges),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'grid3d%bound(ib)%be2v',errorFlag)
    END IF ! errorFlag    

    ALLOCATE(grid3d%bound(ib)%be2c(grid3d%bound(ib)%nBEdges),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'grid3d%bound(ib)%be2c',errorFlag)
    END IF ! errorFlag   
         
    grid3d%bound(ib)%nBEdges = 0
  
    DO ie = 1,grid3d%bound(ib)%nEdges    
      IF ( grid3d%bound(ib)%e2c(1,ie) == 0 .OR. & 
           grid3d%bound(ib)%e2c(2,ie) == 0 ) THEN 
        grid3d%bound(ib)%nBEdges = grid3d%bound(ib)%nBEdges + 1
        
        grid3d%bound(ib)%be2v(1,grid3d%bound(ib)%nBEdges) = & 
          grid3d%bound(ib)%e2v(1,ie)
        grid3d%bound(ib)%be2v(2,grid3d%bound(ib)%nBEdges) = & 
          grid3d%bound(ib)%e2v(2,ie) 
        grid3d%bound(ib)%be2c(grid3d%bound(ib)%nBEdges) = & 
          MAX(grid3d%bound(ib)%e2c(1,ie),grid3d%bound(ib)%e2c(2,ie))
      END IF ! c1
    END DO ! ie    
  END DO ! ib  
   
! ******************************************************************************
! Compute boundary normals to guide user about which surfaces can be selected.
! At present, only allow z-surfaces to be selected because extrusion is limited
! to z-direction. 
! ******************************************************************************

  WRITE(STDOUT,'(3X,A)') & 
    'Determining planarity and orientation of boundaries...'
!  WRITE(STDOUT,'(3X,A)') 'Enter tolerance:'
!  READ(STDIN,*) toler
  toler = 1.0D-10   
  WRITE(STDOUT,'(3X,A,1X,E13.6)') 'Tolerance:',toler  

  WRITE(STDOUT,'(6X,A,5X,A,5(9X,A),7X,A,2(1X,A),3X,A)') '#','nxmin','nxmax', &
        'nymin','nymax','nzmin','nzmax','xflat','yflat','zflat','normDir'
  
  DO ib = 1,grid3d%nBounds
    nxmin =  HUGE(1.0D0)
    nymin =  HUGE(1.0D0)
    nzmin =  HUGE(1.0D0)
    nxmax = -HUGE(1.0D0)
    nymax = -HUGE(1.0D0)
    nzmax = -HUGE(1.0D0)              
  
    grid3d%bound(ib)%selectFlag = .FALSE.

! ==============================================================================
!   Triangles
! ==============================================================================  
  
    DO it = 1,grid3d%bound(ib)%nTris
      v1 = grid3d%bound(ib)%tri2v(1,it)
      v2 = grid3d%bound(ib)%tri2v(2,it)
      v3 = grid3d%bound(ib)%tri2v(3,it)            
      
      v1g = grid3d%bound(ib)%v(v1)
      v2g = grid3d%bound(ib)%v(v2)
      v3g = grid3d%bound(ib)%v(v3)
      
      x1 = gridCENTAUR%xyz(1,v1g)
      y1 = gridCENTAUR%xyz(2,v1g)
      z1 = gridCENTAUR%xyz(3,v1g)
      x2 = gridCENTAUR%xyz(1,v2g)
      y2 = gridCENTAUR%xyz(2,v2g)
      z2 = gridCENTAUR%xyz(3,v2g)
      x3 = gridCENTAUR%xyz(1,v3g)
      y3 = gridCENTAUR%xyz(2,v3g)
      z3 = gridCENTAUR%xyz(3,v3g) 
            
      nx =  0.5D0*((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))
      ny = -0.5D0*((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2))
      nz =  0.5D0*((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2))
      
      nmag = SQRT(nx*nx + ny*ny + nz*nz)
      
      nx = nx/nmag
      ny = ny/nmag
      nz = nz/nmag
      
      nxmin = MIN(nx,nxmin)
      nymin = MIN(ny,nymin)
      nzmin = MIN(nz,nzmin) 
      
      nxmax = MAX(nx,nxmax)
      nymax = MAX(ny,nymax)
      nzmax = MAX(nz,nzmax)                
    END DO ! it

! ==============================================================================
!   Triangles
! ==============================================================================  
  
    DO iq = 1,grid3d%bound(ib)%nQuads
      v1 = grid3d%bound(ib)%quad2v(1,iq)
      v2 = grid3d%bound(ib)%quad2v(2,iq)
      v3 = grid3d%bound(ib)%quad2v(3,iq)
      v4 = grid3d%bound(ib)%quad2v(4,iq)                  
      
      v1g = grid3d%bound(ib)%v(v1)
      v2g = grid3d%bound(ib)%v(v2)
      v3g = grid3d%bound(ib)%v(v3)
      v4g = grid3d%bound(ib)%v(v4)      
      
      x1 = gridCENTAUR%xyz(1,v1g)
      y1 = gridCENTAUR%xyz(2,v1g)
      z1 = gridCENTAUR%xyz(3,v1g)
      x2 = gridCENTAUR%xyz(1,v2g)
      y2 = gridCENTAUR%xyz(2,v2g)
      z2 = gridCENTAUR%xyz(3,v2g)
      x3 = gridCENTAUR%xyz(1,v3g)
      y3 = gridCENTAUR%xyz(2,v3g)
      z3 = gridCENTAUR%xyz(3,v3g) 
      x4 = gridCENTAUR%xyz(1,v4g)
      y4 = gridCENTAUR%xyz(2,v4g)
      z4 = gridCENTAUR%xyz(3,v4g) 

      nx = 0.5D0*((z3 - z1)*(y2 - y4)-(y3 - y1)*(z2 - z4))
      ny = 0.5D0*((x3 - x1)*(z2 - z4)-(z3 - z1)*(x2 - x4))
      nz = 0.5D0*((y3 - y1)*(x2 - x4)-(x3 - x1)*(y2 - y4))

      nmag = SQRT(nx*nx + ny*ny + nz*nz)
      
      nx = nx/nmag
      ny = ny/nmag
      nz = nz/nmag
                        
      nxmin = MIN(nx,nxmin)
      nymin = MIN(ny,nymin)
      nzmin = MIN(nz,nzmin) 
      
      nxmax = MAX(nx,nxmax)
      nymax = MAX(ny,nymax)
      nzmax = MAX(nz,nzmax)                
    END DO ! iq  
  
! ==============================================================================
!   Check for planarity
! ==============================================================================  
    
    IF ( SIGN(1.0D0,nxmin) == SIGN(1.0D0,nxmax) ) THEN 
      IF ( ABS(nxmax-nxmin) < toler ) THEN 
        xFlatFlag = .TRUE.   
      ELSE 
        xFlatFlag = .FALSE. 
      END IF ! ABS
    ELSE 
      IF ( (ABS(nxmin) < toler) .AND. & 
           (ABS(nxmax) < toler) ) THEN 
        xFlatFlag = .TRUE.            
      ELSE 
        xFlatFlag = .FALSE. 
      END IF ! ABS
    END IF ! SIGN
    
    IF ( SIGN(1.0D0,nymin) == SIGN(1.0D0,nymax) ) THEN 
      IF ( ABS(nymax-nymin) < toler ) THEN 
        yFlatFlag = .TRUE.   
      ELSE 
        yFlatFlag = .FALSE. 
      END IF ! ABS
    ELSE 
      IF ( (ABS(nymin) < toler) .AND. & 
           (ABS(nymax) < toler) ) THEN 
        yFlatFlag = .TRUE.            
      ELSE 
        yFlatFlag = .FALSE. 
      END IF ! ABS
    END IF ! SIGN      
    
    IF ( SIGN(1.0D0,nzmin) == SIGN(1.0D0,nzmax) ) THEN 
      IF ( ABS(nzmax-nzmin) < toler ) THEN 
        zFlatFlag = .TRUE.   
      ELSE 
        zFlatFlag = .FALSE. 
      END IF ! ABS
    ELSE 
      IF ( (ABS(nzmin) < toler) .AND. & 
           (ABS(nzmax) < toler) ) THEN 
        zFlatFlag = .TRUE.            
      ELSE 
        zFlatFlag = .FALSE. 
      END IF ! ABS
    END IF ! SIGN  

! ==============================================================================
!   Check for alignment with coordinate axes
! ==============================================================================  
        
    surfNormDir = 0   
    
    IF ( (xFlatFlag .EQV. .TRUE.) .AND. & 
         (yFlatFlag .EQV. .TRUE.) .AND. &
         (zFlatFlag .EQV. .TRUE.) ) THEN 
      term = MAX(ABS(nxmin),ABS(nxmax), & 
                 ABS(nymin),ABS(nymax), &
                 ABS(nzmin),ABS(nzmax))                  
      IF ( (ABS(term) - 1.0D0) < toler ) THEN 
        IF ( term == ABS(nxmin) .OR. term == ABS(nxmax) ) THEN 
          surfNormDir = 1
        ELSE IF ( term == ABS(nymin) .OR. term == ABS(nymax) ) THEN
          surfNormDir = 2        
        ELSE IF ( term == ABS(nzmin) .OR. term == ABS(nzmax) ) THEN
          surfNormDir = 3       
        END IF ! term
      END IF ! ABS            
    END IF ! xFlatFlag

    WRITE(STDOUT,'(5X,I2,6(1X,E13.6),3(1X,L1),1X,I1)') & 
          ib,nxmin,nxmax,nymin,nymax,nzmin,nzmax, & 
          xFlatFlag,yFlatFlag,zFlatFlag,surfNormDir

    IF ( surfNormDir == 3 .AND. nzmax == 1.0D0 ) THEN 
      grid3d%bound(ib)%selectFlag = .TRUE. 
    END IF ! surfNormDir       
  END DO ! ib  

! ******************************************************************************
! Count number of boundaries which can be selected. Store latest one in ib2d,
! which means that if there is only one boundary, have that boundary selected 
! automatically.
! ******************************************************************************

  selectFlagCntr = 0

  DO ib = 1,grid3d%nBounds
    IF ( grid3d%bound(ib)%selectFlag .EQV. .TRUE. ) THEN 
      selectFlagCntr = selectFlagCntr + 1
      ib2d = ib
    END IF ! 
  END DO ! ib

  WRITE(STDOUT,'(5X,A,1X,I2,1X,A)') 'Detected',selectFlagCntr, &
    'boundary/boundaries which can be used for extrusion.' 

! ******************************************************************************
! Select boundary which is to become surface to be extruded
! ******************************************************************************

  IF ( selectFlagCntr == 1 ) THEN 
    WRITE(STDOUT,'(5X,A)') & 
      'Only 1 boundary can be used for extrusion, will be chosen automatically.'
    WRITE(STDOUT,'(5X,A,1X,I2,1X,A)') & 
      'Boundary',ib2d,'will be chosen for extrusion.'             
  ELSE IF ( selectFlagCntr > 1 ) THEN 
    WRITE(STDOUT,'(5X,A)') 'Multiple boundaries can be used for extrusion.' 
    
    DO ib = 1,gridCENTAUR%nBounds
      IF ( ib > 1 ) THEN 
        nBTris  = gridCENTAUR%bInfo(2,ib) - gridCENTAUR%bInfo(2,ib-1) 
        nBQuads = gridCENTAUR%bInfo(3,ib) - gridCENTAUR%bInfo(3,ib-1) 
      ELSE 
        nBTris  = gridCENTAUR%bInfo(2,ib)
        nBQuads = gridCENTAUR%bInfo(3,ib)
      END IF ! ib 

      WRITE(STDOUT,'(5X,I2,1X,A20,1X,I4,2(1X,I5),1X,L1)') & 
        ib,ADJUSTL(TRIM(gridCENTAUR%bName(ib))),gridCENTAUR%bInfo(1,ib), & 
        nBTris,nBQuads,grid3d%bound(ib)%selectFlag
    END DO ! ib
        
    choiceLoop: DO               
      WRITE(STDOUT,'(3X,A)') 'Enter boundary which becomes 2d grid:'
      READ(STDIN,*) ib2d  

      IF ( grid3d%bound(ib2d)%selectFlag .EQV. .TRUE. ) THEN
        EXIT choiceLoop
      ELSE 
        WRITE(STDOUT,'(3X,A)') 'This boundary may not be selected! Try again...'
      END IF ! grid3d%bound(ib2d)%selectFlag
    END DO choiceLoop
  ELSE 
    WRITE(STDOUT,'(5X,A)') 'It appears no boundaries can be used for extrusion.'
    WRITE(STDOUT,'(5X,A)') 'Select manually if you know what you are doing...'
    WRITE(STDOUT,'(3X,A)') 'Enter boundary which becomes 2d grid:'
    READ(STDIN,*) ib2d                    
  END IF ! selectFlagCntr 
   
! ******************************************************************************
! Identify shared boundaries and edges
! ******************************************************************************

  WRITE(STDOUT,'(3X,A)') 'Identifying shared boundary edges...'

  ALLOCATE(bCntr(grid3d%nBounds ),STAT=errorFlag)      
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'bCntr',errorFlag)
  END IF ! errorFlag             

  DO ib = 1,grid3d%nBounds  
    bCntr(ib) = 0
  END DO ! ib  
              
  ALLOCATE(bPntr(2,grid3d%bound(ib2d)%nBEdges),STAT=errorFlag)      
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'bPntr',errorFlag)
  END IF ! errorFlag             
                
  DO ie = 1,grid3d%bound(ib2d)%nBEdges
    bPntr(1,ie) = 0   
    bPntr(2,ie) = 0
  END DO ! ie                
                
  DO ib = 1,grid3d%nBounds  
    IF ( ib /= ib2d ) THEN 
      DO ie = 1,grid3d%bound(ib)%nBEdges
        v1 = grid3d%bound(ib)%be2v(1,ie)
        v2 = grid3d%bound(ib)%be2v(2,ie)
        
        v1g = grid3d%bound(ib)%v(v1)
        v2g = grid3d%bound(ib)%v(v2)                        
                
        DO ie2 = 1,grid3d%bound(ib2d)%nBEdges
          v12 = grid3d%bound(ib2d)%be2v(1,ie2)
          v22 = grid3d%bound(ib2d)%be2v(2,ie2)

          v12g = grid3d%bound(ib2d)%v(v12)
          v22g = grid3d%bound(ib2d)%v(v22)        
                    
          IF ( v1g == v22g .AND. v2g == v12g ) THEN
            bPntr(1,ie2) = ib
            bPntr(2,ie2) = ie            
          END IF ! v1          
        END DO ! ie2
      END DO ! ie    
    END IF ! ib  
  END DO ! ib 
      
  DO ie = 1,grid3d%bound(ib2d)%nBEdges
    bCntr(bPntr(1,ie)) = bCntr(bPntr(1,ie)) + 1
  END DO ! ie         

! ******************************************************************************
! Set up boundary data structure for 2d grid
! ******************************************************************************

! ==============================================================================
! Count number of boundaries
! ==============================================================================

  grid%nBounds = 0

  DO ib = 1,grid3d%nBounds  
    WRITE(STDOUT,'(5X,I2,1X,I5)') ib,bCntr(ib)
    
    IF ( bCntr(ib) /= 0 ) THEN 
      grid%nBounds = grid%nBounds + 1
    END IF ! bCntr
  END DO ! ib  

  WRITE(STDOUT,'(3X,A,1X,I2,1X,A)') 'Detected',grid%nBounds,'boundaries.'
              
  ALLOCATE(grid%bound(grid%nBounds),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound',errorFlag)
  END IF ! errorFlag

              
  ALLOCATE(grid3d%boundMap(grid3d%nBounds),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound',errorFlag)
  END IF ! errorFlag

! ==============================================================================
! Set up mapping from 3d to 2d boundaries
! ==============================================================================

  WRITE(STDOUT,'(3X,A)') 'Generating boundary mapping...'

  ib = 0

  DO ib3d = 1,grid3d%nBounds      
    IF ( bCntr(ib3d) /= 0 ) THEN
      ib = ib + 1 
      grid3d%boundMap(ib3d) = ib 
    ELSE 
      grid3d%boundMap(ib3d) = 0     
    END IF ! bCntr
    
    WRITE(STDOUT,'(5X,I2,1X,I5)') ib3d,grid3d%boundMap(ib3d)
  END DO ! ib3d

! ==============================================================================
! Specify number of edges and boundary type for 2d grid
! ==============================================================================

  ib = 0

  DO ib3d = 1,grid3d%nBounds      
    IF ( grid3d%boundMap(ib3d) /= 0 ) THEN
      ib = grid3d%boundMap(ib3d)
    
      grid%bound(ib)%nEdges = bCntr(ib3d)
      grid%bound(ib)%bName  = gridCENTAUR%bName(ib3d)
      grid%bound(ib)%bType  = gridCENTAUR%bInfo(1,ib3d)
            
      ALLOCATE(grid%bound(ib)%e2v(2,grid%bound(ib)%nEdges),STAT=errorFlag)
      IF ( errorFlag /= NO_ERROR ) THEN
        CALL errorHandling(ALLOCATE_ERROR,'grid%bound%e2v',errorFlag)
      END IF ! errorFlag      
    END IF ! grid3d%boundMap(ib3d)
  END DO ! ib

! ==============================================================================
! Specify edge connectivity for 2d grid
! ==============================================================================

  DO ib = 1,grid%nBounds    
    ALLOCATE(grid%bound(ib)%e2v(2,grid%bound(ib)%nEdges),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN
      CALL errorHandling(ALLOCATE_ERROR,'grid%bound%e2v',errorFlag)
    END IF ! errorFlag
  END DO ! ib

  DO ib = 1,grid%nBounds
    grid%bound(ib)%nEdges = 0 ! Reset
  END DO ! ib

  DO ie2 = 1,grid3d%bound(ib2d)%nBEdges
    ib3d = bPntr(1,ie2)
    
    IF ( grid3d%boundMap(ib3d) /= 0 ) THEN 
      ib = grid3d%boundMap(ib3d)
      
      grid%bound(ib)%nEdges = grid%bound(ib)%nEdges + 1
        
      grid%bound(ib)%e2v(1,grid%bound(ib)%nEdges) = & 
        grid3d%bound(ib2d)%be2v(1,ie2)
      grid%bound(ib)%e2v(2,grid%bound(ib)%nEdges) = & 
        grid3d%bound(ib2d)%be2v(2,ie2)
    END IF ! grid3d%boundMap(ib3d)            
  END DO ! ie

! ==============================================================================
! Store 2d connectivity
! ==============================================================================

  grid%nTris  = grid3d%bound(ib2d)%nTris 
  grid%nQuads = grid3d%bound(ib2d)%nQuads
  grid%nCells = grid%nTris + grid%nQuads
  
  IF ( grid%nTris > 0 ) THEN 
    ALLOCATE(grid%tri2v(3,grid%nTris),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'grid%tri2v',errorFlag)
    END IF ! errorFlag

    DO it = 1,grid3d%bound(ib2d)%nTris
      grid%tri2v(1,it) = grid3d%bound(ib2d)%tri2v(1,it)
      grid%tri2v(2,it) = grid3d%bound(ib2d)%tri2v(2,it)
      grid%tri2v(3,it) = grid3d%bound(ib2d)%tri2v(3,it)       
    END DO ! it
  END IF ! grid%nTris  
   
  IF ( grid%nQuads > 0 ) THEN 
    ALLOCATE(grid%quad2v(4,grid%nQuads),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'grid%quad2v',errorFlag)
    END IF ! errorFlag

    DO iq = 1,grid3d%bound(ib2d)%nQuads
      grid%quad2v(1,iq) = grid3d%bound(ib2d)%quad2v(1,iq)
      grid%quad2v(2,iq) = grid3d%bound(ib2d)%quad2v(2,iq)
      grid%quad2v(3,iq) = grid3d%bound(ib2d)%quad2v(3,iq)
      grid%quad2v(4,iq) = grid3d%bound(ib2d)%quad2v(4,iq)             
    END DO ! iq
  END IF ! grid%nTris     

! ==============================================================================
! Store 2d coordinates
! ==============================================================================

  grid%nVert = grid3d%bound(ib2d)%nVert   
   
  ALLOCATE(grid%xy(2,grid%nVert),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'grid%xy',errorFlag)
  END IF ! errorFlag   
   
  DO iv = 1,grid%nVert
    ivg = grid3d%bound(ib2d)%v(iv)
    
    grid%xy(1,iv) = gridCENTAUR%xyz(1,ivg)
    grid%xy(2,iv) = gridCENTAUR%xyz(2,ivg)      
  END DO ! iv 
   
! ******************************************************************************
! Destroy CENTAUR data because will need to be reused for converting back to 3d
! ******************************************************************************
   
  IF ( gridCENTAUR%nTets > 0 ) THEN 
    DEALLOCATE(gridCENTAUR%tet2v,STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(DEALLOCATE_ERROR,'gridCENTAUR%tet2v',errorFlag)
    END IF ! errorFlag
  END IF ! gridCENTAUR   
   
  IF ( gridCENTAUR%nHexs > 0 ) THEN 
    DEALLOCATE(gridCENTAUR%hex2v,STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(DEALLOCATE_ERROR,'gridCENTAUR%hex2v',errorFlag)
    END IF ! errorFlag
  END IF ! gridCENTAUR

  IF ( gridCENTAUR%nPris > 0 ) THEN 
    DEALLOCATE(gridCENTAUR%pri2v,STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(DEALLOCATE_ERROR,'gridCENTAUR%pri2v',errorFlag)
    END IF ! errorFlag
  END IF ! gridCENTAUR

  IF ( gridCENTAUR%nPyrs > 0 ) THEN 
    DEALLOCATE(gridCENTAUR%pyr2v,STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(DEALLOCATE_ERROR,'gridCENTAUR%pyr2v',errorFlag)
    END IF ! errorFlag
  END IF ! gridCENTAUR  
  
  DEALLOCATE(gridCENTAUR%xyz,STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(DEALLOCATE_ERROR,'gridCENTAUR%xyz',errorFlag)
  END IF ! errorFlag

  IF ( gridCENTAUR%nBTris > 0 ) THEN 
    DEALLOCATE(gridCENTAUR%bTri2v,STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(DEALLOCATE_ERROR,'gridCENTAUR%bTri2v',errorFlag)
    END IF ! errorFlag   
  END IF ! gridCENTAUR

  IF ( gridCENTAUR%nBQuads > 0 ) THEN
    DEALLOCATE(gridCENTAUR%bQuad2v,STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(DEALLOCATE_ERROR,'gridCENTAUR%bQuad2v',errorFlag)
    END IF ! errorFlag  
  END IF ! gridCENTAUR

  DEALLOCATE(gridCENTAUR%bInfo,STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(DEALLOCATE_ERROR,'gridCENTAUR%bInfo',errorFlag)
  END IF ! errorFlag

  DEALLOCATE(gridCENTAUR%bName,STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(DEALLOCATE_ERROR,'gridCENTAUR%bName',errorFlag)
  END IF ! errorFlag
    
  gridCENTAUR%nBounds = 0
  gridCENTAUR%nBTris  = 0
  gridCENTAUR%nBQuads = 0 
  gridCENTAUR%nCells  = 0
  gridCENTAUR%nHexs   = 0
  gridCENTAUR%nPris   = 0
  gridCENTAUR%nTets   = 0
  gridCENTAUR%nPyrs   = 0
  gridCENTAUR%nVert   = 0
   
! ******************************************************************************
! End
! ******************************************************************************

  WRITE(STDOUT,'(1X,A)') 'Creating 2d-grid in generic format done.'
  
END SUBROUTINE CENTAUR2Generic
