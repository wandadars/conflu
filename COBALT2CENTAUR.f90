! ******************************************************************************
!
! $Id: COBALT2CENTAUR.f90,v 1.1 2005/03/05 17:25:53 haselbac Exp $
!
! Filename: COBALT2CENTAUR.F90
!
! Purpose: Convert COBALT grid format into CENTAUR grid format.
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
! Copyright: (c) 2002 by the University of Illinois
!
! RCS Revision history:
!
!   $Log: COBALT2CENTAUR.f90,v $
!   Revision 1.1  2005/03/05 17:25:53  haselbac
!   Initial revision
!
!   Revision 1.1  2003/03/07 15:26:34  haselbac
!   Initial revision
!
!
! ******************************************************************************

SUBROUTINE COBALT2CENTAUR

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

  INTEGER :: c1,c2,i,iBoundIndex,ic,j,k,l,flag,nBoundsCounter,v1,v2,v3,v4
  INTEGER, DIMENSION(3) :: pyr2vTemp
  INTEGER, DIMENSION(:), ALLOCATABLE :: cellType,cellMapp
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: cellDeg,nBFaces
  DOUBLE PRECISION :: cpx,cpy,cpz,dx,dy,dz,xm,x1,x2,x3,x4,ym,y1,y2,y3,y4, & 
                      zm,z1,z2,z3,z4
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: cofg

! ******************************************************************************
! Start
! ******************************************************************************

  WRITE(STDOUT,'(/,1X,A)') 'Converting from COBALT to CENTAUR format...'
   
  
! ******************************************************************************
! Determine number of different cell types 
! ******************************************************************************

  WRITE(STDOUT,'(3X,A)') 'Determining number of different cell types...'

  ALLOCATE(cellDeg(FACE_TYPE_TRI:FACE_TYPE_QUAD,gridCENTAUR%nCells),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'cellDeg',errorFlag)
  END IF ! errorFlag

  cellDeg(:,:) = 0

  DO i = 1,gridCOBALT%nFaces
    c1 = gridCOBALT%f2c(1,i)
    c2 = gridCOBALT%f2c(2,i)

    IF ( gridCOBALT%nvpf(i) == 3 ) THEN ! triangular face 
      IF ( c1 > 0 ) THEN 
        cellDeg(FACE_TYPE_TRI,c1) = cellDeg(FACE_TYPE_TRI,c1) + 1
      END IF ! c1

      IF ( c2 > 0 ) THEN 
        cellDeg(FACE_TYPE_TRI,c2) = cellDeg(FACE_TYPE_TRI,c2) + 1
      END IF ! c2
    ELSE IF ( gridCOBALT%nvpf(i) == 4 ) THEN ! quadrilateral face
      IF ( c1 > 0 ) THEN
        cellDeg(FACE_TYPE_QUAD,c1) = cellDeg(FACE_TYPE_QUAD,c1) + 1
      END IF ! c1

      IF ( c2 > 0 ) THEN
        cellDeg(FACE_TYPE_QUAD,c2) = cellDeg(FACE_TYPE_QUAD,c2) + 1
      END IF ! c2
    ELSE ! should be impossible due to previous check
      CALL errorHandling(INVALID_FACETYPE_ERROR)
    END IF ! gridCOBALT    
  END DO ! i

  ALLOCATE(cellType(gridCENTAUR%nCells),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'cellType',errorFlag)
  END IF ! errorFlag

  ALLOCATE(cellMapp(gridCENTAUR%nCells),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'cellMapp',errorFlag)
  END IF ! errorFlag

  cellType(:) = 0
  cellMapp(:) = 0

  gridCENTAUR%nTets = 0
  gridCENTAUR%nHexs = 0
  gridCENTAUR%nPris = 0
  gridCENTAUR%nPyrs = 0

  DO i = 1,gridCENTAUR%nCells
    IF ( cellDeg(FACE_TYPE_QUAD,i) == 0 .AND. & 
         cellDeg(FACE_TYPE_TRI, i) == 4 ) THEN ! Tetrahedron
      gridCENTAUR%nTets = gridCENTAUR%nTets + 1
      cellType(i) = CELL_TYPE_TETRAHEDRON
      cellMapp(i) = gridCENTAUR%nTets
    ELSE IF ( cellDeg(FACE_TYPE_QUAD,i) == 1 .AND. & 
              cellDeg(FACE_TYPE_TRI, i) == 4 ) THEN ! Pyramid
      gridCENTAUR%nPyrs = gridCENTAUR%nPyrs + 1
      cellType(i) = CELL_TYPE_PYRAMID 
      cellMapp(i) = gridCENTAUR%nPyrs         
    ELSE IF ( cellDeg(FACE_TYPE_QUAD,i) == 3 .AND. & 
              cellDeg(FACE_TYPE_TRI, i) == 2 ) THEN ! Prism
      gridCENTAUR%nPris = gridCENTAUR%nPris + 1
      cellType(i) = CELL_TYPE_PRISM 
      cellMapp(i) = gridCENTAUR%nPris           
    ELSE IF ( cellDeg(FACE_TYPE_QUAD,i) == 6 .AND. & 
              cellDeg(FACE_TYPE_TRI, i) == 0 ) THEN ! Hexahedron
      gridCENTAUR%nHexs = gridCENTAUR%nHexs + 1
      cellType(i) = CELL_TYPE_HEXAHEDRON
      cellMapp(i) = gridCENTAUR%nHexs      
    ELSE 
      CALL errorHandling(INVALID_CELLTYPE_ERROR)
    END IF ! cellDeg
  END DO ! i

  IF ( gridCENTAUR%nTets > 0 ) THEN 
    ALLOCATE(gridCENTAUR%tet2v(4,gridCENTAUR%nTets),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%tet2v',errorFlag)
    END IF ! errorFlag 
  END IF ! gridCENTAUR%nTets

  IF ( gridCENTAUR%nHexs > 0 ) THEN 
    ALLOCATE(gridCENTAUR%hex2v(8,gridCENTAUR%nHexs),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%hex2v',errorFlag)
    END IF ! errorFlag      
  END IF ! gridCENTAUR%nHexs

  IF ( gridCENTAUR%nPris > 0 ) THEN 
    ALLOCATE(gridCENTAUR%pri2v(6,gridCENTAUR%nPris),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%pri2v',errorFlag)
    END IF ! errorFlag  
  END IF ! gridCENTAUR%nPris

  IF ( gridCENTAUR%nPyrs > 0 ) THEN 
    ALLOCATE(gridCENTAUR%pyr2v(5,gridCENTAUR%nPyrs),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%pyr2v',errorFlag)
    END IF ! errorFlag
  END IF ! gridCENTAUR%nPyrs

  WRITE(STDOUT,'(3X,A)') 'Cell counts:'
  WRITE(STDOUT,'(5X,A,1X,I9)') 'Tetrahedra:',gridCENTAUR%nTets
  WRITE(STDOUT,'(5X,A,1X,I9)') 'Hexahedra: ',gridCENTAUR%nHexs
  WRITE(STDOUT,'(5X,A,1X,I9)') 'Prisms:    ',gridCENTAUR%nPris
  WRITE(STDOUT,'(5X,A,1X,I9)') 'Pyramids:  ',gridCENTAUR%nPyrs            

! ******************************************************************************
! Loop over faces and sort into respective cell types
! ******************************************************************************

  WRITE(STDOUT,'(3X,A)') 'Building cell connectivity...'

  cellDeg(:,:) = 0

  DO i = 1,gridCOBALT%nFaces
    DO j = 1,2
      c1 = gridCOBALT%f2c(j,i)
    
! ==============================================================================
!     Triangular faces
! ==============================================================================    
    
!      IF ( gridCOBALT%nvpf(i) == 3 ) THEN ! triangular face
        IF ( c1 > 0 ) THEN ! cell interior to solution domain
          SELECT CASE ( cellType(c1) ) 
          
! ------------------------------------------------------------------------------
!           Tetrahedra
! ------------------------------------------------------------------------------         
          
            CASE ( CELL_TYPE_TETRAHEDRON )     
              IF ( cellDeg(FACE_TYPE_TRI,c1) < 2 ) THEN ! only need two faces
                cellDeg(FACE_TYPE_TRI,c1) = cellDeg(FACE_TYPE_TRI,c1) + 1

                ic = cellMapp(c1)

                IF ( cellDeg(FACE_TYPE_TRI,c1) == 1 ) THEN 
                  gridCENTAUR%tet2v(1:3,ic) = gridCOBALT%f2v(1:3,i)
                ELSE 
                  DO k = 1,3
                    flag = 0
                    DO l = 1,3
                      IF ( gridCOBALT%f2v(k,i) == gridCENTAUR%tet2v(l,ic) ) THEN 
                        flag = 1
                        EXIT
                      END IF ! gridCOBALT 
                    END DO ! l
                    
                    IF ( flag == 0 ) THEN 
                      EXIT
                    END IF ! flag
                  END DO ! k
                  
                  IF ( flag == 0 ) THEN 
                    gridCENTAUR%tet2v(4,ic) = gridCOBALT%f2v(k,i)
                  END IF ! flag
                END IF ! cellDeg
              END IF ! cellDeg

! ------------------------------------------------------------------------------
!           Prism
! ------------------------------------------------------------------------------         

            CASE ( CELL_TYPE_PRISM ) 
              IF ( gridCOBALT%nvpf(i) == 3 .AND. & 
                   cellDeg(FACE_TYPE_TRI,c1) < 2 ) THEN ! only need 2 tri faces
                cellDeg(FACE_TYPE_TRI,c1) = cellDeg(FACE_TYPE_TRI,c1) + 1             

                ic = cellMapp(c1)

                IF ( cellDeg(FACE_TYPE_TRI,c1) == 1 ) THEN 
                  gridCENTAUR%pri2v(1:3,ic) = gridCOBALT%f2v(1:3,i)
                ELSE IF ( cellDeg(FACE_TYPE_TRI,c1) == 2 ) THEN 
                  gridCENTAUR%pri2v(4:6,ic) = gridCOBALT%f2v(1:3,i)
                ELSE  
                  WRITE(*,*) 'ERROR - more than 2 tri faces in prism' 
                  STOP
                END IF ! cellDeg                
              END IF ! cellDeg

! ------------------------------------------------------------------------------
!           Pyramid
! ------------------------------------------------------------------------------         

            CASE ( CELL_TYPE_PYRAMID )      
              IF ( gridCOBALT%nvpf(i) == 3 .AND. & 
                   cellDeg(FACE_TYPE_TRI,c1) < 1 ) THEN ! only need 1 tri face
                cellDeg(FACE_TYPE_TRI,c1) = cellDeg(FACE_TYPE_TRI,c1) + 1

                ic = cellMapp(c1)

                IF ( cellDeg(FACE_TYPE_QUAD,c1) == 0 ) THEN 
                  gridCENTAUR%pyr2v(1:3,ic) = gridCOBALT%f2v(1:3,i)                               
                ELSE 
                  DO k = 1,3
                    flag = 0
                    DO l = 1,4
                      IF ( gridCOBALT%f2v(k,i) == gridCENTAUR%pyr2v(l,ic) ) THEN 
                        flag = 1
                        EXIT
                      END IF ! gridCOBALT 
                    END DO ! l
                           
                    IF ( flag == 0 ) THEN
                      EXIT
                    END IF ! flag
                  END DO ! k
                  
                  IF ( flag == 0 ) THEN 
                    gridCENTAUR%pyr2v(5,ic) = gridCOBALT%f2v(k,i)
                  END IF ! flag
                END IF ! cellDeg
              END IF ! cellDeg

              IF ( gridCOBALT%nvpf(i) == 4 .AND. & 
                   cellDeg(FACE_TYPE_QUAD,c1) < 1 ) THEN ! only need one face
                cellDeg(FACE_TYPE_QUAD,c1) = cellDeg(FACE_TYPE_QUAD,c1) + 1

                ic = cellMapp(c1)

                IF ( cellDeg(FACE_TYPE_TRI,c1) == 0 ) THEN 
                  gridCENTAUR%pyr2v(1:4,ic) = gridCOBALT%f2v(1:4,i)                       
                ELSE 
                  pyr2vTemp(1:3) = gridCENTAUR%pyr2v(1:3,ic)
                  gridCENTAUR%pyr2v(1:4,ic) = gridCOBALT%f2v(1:4,i)
                
                  DO k = 1,3
                    flag = 0
                    DO l = 1,4
                      IF ( pyr2vTemp(k) == gridCENTAUR%pyr2v(l,ic) ) THEN 
                        flag = 1
                        EXIT
                      END IF ! gridCOBALT 
                    END DO ! l
                     
                    IF ( flag == 0 ) THEN
                      EXIT
                    END IF ! flag
                  END DO ! k
                  
                  IF ( flag == 0 ) THEN 
                    gridCENTAUR%pyr2v(5,ic) = pyr2vTemp(k)                  
                  END IF ! flag
                END IF ! cellDeg
              END IF ! cellDeg

! ------------------------------------------------------------------------------
!           Default
! ------------------------------------------------------------------------------         

            CASE DEFAULT ! hexahedral cell impossible
              CALL errorHandling(REACHED_DEFAULT)
          END SELECT ! cellType
        ELSE ! cell exterior to solution domain
      
        END IF ! c1
        
! ==============================================================================
!     Quadrilateral faces
! ==============================================================================    
        
!      ELSE IF ( gridCOBALT%nvpf(i) == 4 ) THEN ! quadrilateral face
!      
!! BEGIN TEMPORARY
!        WRITE(STDOUT,*) 'Reached quad face part of if block.'
!        STOP
!! END TEMPORARY   
!   
!      ELSE ! should be impossible, as already checked - see above
!        CALL errorHandling(INVALID_FACETYPE_ERROR)
!      END IF ! gridCOBALT%nvpf    
    END DO ! j
  END DO ! i

! ******************************************************************************
! Check orientation
! ******************************************************************************

  WRITE(STDOUT,'(3X,A)') 'Checking cell connectivity...'

  ALLOCATE(cofg(3,gridCENTAUR%nCells),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'cofg',errorFlag)
  END IF ! errorFlag  

! ------------------------------------------------------------------------------
! Tetrahedra
! ------------------------------------------------------------------------------

  DO i = 1,gridCENTAUR%nTets
    v1 = gridCENTAUR%tet2v(1,i)
    v2 = gridCENTAUR%tet2v(2,i)
    v3 = gridCENTAUR%tet2v(3,i)
    v4 = gridCENTAUR%tet2v(4,i)            
  
    x1 = gridCENTAUR%xyz(1,v1)
    y1 = gridCENTAUR%xyz(2,v1)
    z1 = gridCENTAUR%xyz(3,v1)        
    x2 = gridCENTAUR%xyz(1,v2)
    y2 = gridCENTAUR%xyz(2,v2)
    z2 = gridCENTAUR%xyz(3,v2) 
    x3 = gridCENTAUR%xyz(1,v3)
    y3 = gridCENTAUR%xyz(2,v3)
    z3 = gridCENTAUR%xyz(3,v3) 
    x4 = gridCENTAUR%xyz(1,v4)
    y4 = gridCENTAUR%xyz(2,v4)
    z4 = gridCENTAUR%xyz(3,v4)         
  
    cofg(1,i) = 0.25D0*(x1 + x2 + x3 + x4)
    cofg(2,i) = 0.25D0*(y1 + y2 + y3 + y4)
    cofg(3,i) = 0.25D0*(z1 + z2 + z3 + z4)        
  
    flag = 0
  
    DO j = 1,4 ! loop over faces
      v1 = gridCENTAUR%tet2v(f2vTet(1,j),i)
      v2 = gridCENTAUR%tet2v(f2vTet(2,j),i)
      v3 = gridCENTAUR%tet2v(f2vTet(3,j),i) 
      
      x1 = gridCENTAUR%xyz(1,v1)
      y1 = gridCENTAUR%xyz(2,v1)
      z1 = gridCENTAUR%xyz(3,v1)        
      x2 = gridCENTAUR%xyz(1,v2)
      y2 = gridCENTAUR%xyz(2,v2)
      z2 = gridCENTAUR%xyz(3,v2) 
      x3 = gridCENTAUR%xyz(1,v3)
      y3 = gridCENTAUR%xyz(2,v3)
      z3 = gridCENTAUR%xyz(3,v3)        
      
      xm = (x1 + x2 + x3)/3.0D0
      ym = (y1 + y2 + y3)/3.0D0
      zm = (z1 + z2 + z3)/3.0D0
      
      dx = xm - cofg(1,i)
      dy = ym - cofg(2,i)
      dz = zm - cofg(3,i)
      
      cpx =   (y2 - y1)*(z3 - z1) - (z2 - z1)*(y3 - y1)
      cpy = -((x2 - x1)*(z3 - z1) - (z2 - z1)*(x3 - x1))
      cpz =   (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1) 
      
      IF ( cpx*dx + cpy*dy + cpz*dz < 0.0D0 ) THEN 
        flag = flag + 1
        
        IF ( flag > 1 ) THEN 
          WRITE(*,*) 'Error in checking of orientation!'
          STOP
        END IF ! flag
        
        gridCENTAUR%tet2v(f2vTet(1,j),i) = v1
        gridCENTAUR%tet2v(f2vTet(2,j),i) = v3
        gridCENTAUR%tet2v(f2vTet(3,j),i) = v2
      END IF ! cpx           
    END DO ! j    
  END DO ! i

! ------------------------------------------------------------------------------
! Prisms
! ------------------------------------------------------------------------------

!  DO i = 1,gridCENTAUR%nPris
!   
!  END DO ! i  

! ------------------------------------------------------------------------------
! Pyramids
! ------------------------------------------------------------------------------

!  DO i = 1,gridCENTAUR%nPyrs
!   
!  END DO ! i  

! ******************************************************************************
! Generate boundary data structure 
! ******************************************************************************

  WRITE(STDOUT,'(3X,A)') 'Building boundary data structure...'

  nBoundsCounter      = 0
  gridCENTAUR%nBTris  = 0
  gridCENTAUR%nBQuads = 0

  ALLOCATE(gridCENTAUR%bInfo(3,gridCENTAUR%nBounds),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%bInfo',errorFlag)
  END IF ! errorFlag

  gridCENTAUR%bInfo(:,:) = 0

  ALLOCATE(gridCENTAUR%bName(gridCENTAUR%nBounds),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%bName',errorFlag)
  END IF ! errorFlag

  WRITE(STDOUT,'(5X,A)') 'Determining boundary face dimensions...'

  DO i = 1,gridCOBALT%nFaces
    DO j = 1,2
      c1 = gridCOBALT%f2c(j,i)
    
      IF ( c1 < 0 ) THEN              
        IF ( nBoundsCounter /= 0 ) THEN
          flag = 0 
          DO k = 1,gridCENTAUR%nBounds    
            IF ( gridCENTAUR%bInfo(1,k) == c1 ) THEN
              flag = 1
              EXIT
            END IF ! bound
          END DO ! k
          
          IF ( flag == 1 ) THEN 
            iBoundIndex = k         
          ELSE 
            nBoundsCounter = nBoundsCounter + 1
            iBoundIndex    = nBoundsCounter
          END IF ! flag
        ELSE 
          nBoundsCounter = 1
          iBoundIndex    = 1      
        END IF ! nBoundsCounter
                                
        gridCENTAUR%bInfo(1,iBoundIndex) = c1
                                
        IF ( gridCOBALT%nvpf(i) == 3 ) THEN 
          gridCENTAUR%nBTris = gridCENTAUR%nBTris + 1 
          gridCENTAUR%bInfo(2,iBoundIndex) = & 
            gridCENTAUR%bInfo(2,iBoundIndex) + 1 
        ELSE 
          gridCENTAUR%nBQuads = gridCENTAUR%nBQuads + 1
          gridCENTAUR%bInfo(3,iBoundIndex) = & 
            gridCENTAUR%bInfo(3,iBoundIndex) + 1          
        END IF ! gridCOBALT
      END IF ! c1
    END DO ! j
  END DO ! i        
        
  DO i = 2,gridCENTAUR%nBounds
    gridCENTAUR%bInfo(2,i) = gridCENTAUR%bInfo(2,i) + gridCENTAUR%bInfo(2,i-1)
    gridCENTAUR%bInfo(3,i) = gridCENTAUR%bInfo(3,i) + gridCENTAUR%bInfo(3,i-1)    
  END DO ! i  
    
  WRITE(STDOUT,'(3X,A)') 'Boundary patch information:'  
    
  DO i = 1,gridCENTAUR%nBounds
    WRITE(STDOUT,'(5X,3(1X,I6))') gridCENTAUR%bInfo(1:3,i)
  END DO ! i    
          
! ==============================================================================
! Construct boundary face arrays
! ==============================================================================

  WRITE(STDOUT,'(3X,A)') 'Building and checking boundary face connectivity...'

  DO i = 1,gridCENTAUR%nBounds  
    IF ( gridCENTAUR%nBTris > 0 ) THEN 
      ALLOCATE(gridCENTAUR%bTri2v(3,gridCENTAUR%nBTris),STAT=errorFlag)
      IF ( errorFlag /= NO_ERROR ) THEN 
        CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%bTri2v',errorFlag)
      END IF ! errorFlag
    END IF ! bound
    
    IF ( gridCENTAUR%nBQuads > 0 ) THEN 
      ALLOCATE(gridCENTAUR%bQuad2v(4,gridCENTAUR%nBQuads),STAT=errorFlag)
      IF ( errorFlag /= NO_ERROR ) THEN 
        CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%bQuad2v',errorFlag)
      END IF ! errorFlag
    END IF ! bound   
  END DO ! i

  ALLOCATE(nBFaces(2,gridCENTAUR%nBounds),STAT=errorFlag) ! Counter
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'nBFaces',errorFlag)
  END IF ! errorFlag

  nBFaces(:,:) = 0

! ----------------------------------------------------------------------------
! Loop over boundary faces
! ----------------------------------------------------------------------------

  DO i = 1,gridCOBALT%nFaces
    DO j = 1,2 ! loop over the two cells straddling the face
      c1 = gridCOBALT%f2c(j,i)
    
      IF ( c1 < 0 ) THEN ! this face is a boundary face
        flag = 0 
        DO k = 1,gridCENTAUR%nBounds
          IF ( gridCENTAUR%bInfo(1,k) == c1 ) THEN
            EXIT
          END IF ! bound
        END DO ! k

        c2 = gridCOBALT%f2c(3-j,i) ! get other cell
        
        IF ( c2 < 0 ) THEN ! other cell cannot also be a boundary cell
          WRITE(*,*) 'Cell connectivity error!',i,j,c1,c2
          STOP
        END IF ! c2
        
! ----- Triangular face -------------------------------------------------------
        
        IF ( gridCOBALT%nvpf(i) == 3 ) THEN 
          nBFaces(1,k) = nBFaces(1,k) + 1
         
          IF ( k == 1 ) THEN 
            l = nBFaces(1,k)
          ELSE 
            l = gridCENTAUR%bInfo(2,k-1) + nBFaces(1,k)
          END IF ! k
          
          gridCENTAUR%bTri2v(1:3,l) = gridCOBALT%f2v(1:3,i)
          
          v1 = gridCENTAUR%bTri2v(1,l)
          v2 = gridCENTAUR%bTri2v(2,l)
          v3 = gridCENTAUR%bTri2v(3,l)
          
          x1 = gridCENTAUR%xyz(1,v1)
          y1 = gridCENTAUR%xyz(2,v1)
          z1 = gridCENTAUR%xyz(3,v1)        
          x2 = gridCENTAUR%xyz(1,v2)
          y2 = gridCENTAUR%xyz(2,v2)
          z2 = gridCENTAUR%xyz(3,v2) 
          x3 = gridCENTAUR%xyz(1,v3)
          y3 = gridCENTAUR%xyz(2,v3)
          z3 = gridCENTAUR%xyz(3,v3)        

          xm = (x1 + x2 + x3)/3.0D0
          ym = (y1 + y2 + y3)/3.0D0
          zm = (z1 + z2 + z3)/3.0D0

          dx = xm - cofg(1,c2)
          dy = ym - cofg(2,c2)
          dz = zm - cofg(3,c2)

          cpx =   (y2 - y1)*(z3 - z1) - (z2 - z1)*(y3 - y1)
          cpy = -((x2 - x1)*(z3 - z1) - (z2 - z1)*(x3 - x1))
          cpz =   (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1) 

          IF ( cpx*dx + cpy*dy + cpz*dz < 0.0D0 ) THEN 
            gridCENTAUR%bTri2v(1,l) = v1
            gridCENTAUR%bTri2v(2,l) = v3
            gridCENTAUR%bTri2v(3,l) = v2
          END IF ! cpx                    
          
! ----- Quadrilateral face ----------------------------------------------------           
          
        ELSE 
          nBFaces(2,k) = nBFaces(2,k) + 1
          
          IF ( k == 1 ) THEN 
            l = nBFaces(2,k)
          ELSE 
            l = gridCENTAUR%bInfo(3,k-1) + nBFaces(2,k)
          END IF ! k

          gridCENTAUR%bQuad2v(1:4,l) = gridCOBALT%f2v(1:4,i)
          
          v1 = gridCENTAUR%bQuad2v(1,l) ! get first triangular subface
          v2 = gridCENTAUR%bQuad2v(2,l)
          v3 = gridCENTAUR%bQuad2v(3,l)                   
          
          x1 = gridCENTAUR%xyz(1,v1)
          y1 = gridCENTAUR%xyz(2,v1)
          z1 = gridCENTAUR%xyz(3,v1)        
          x2 = gridCENTAUR%xyz(1,v2)
          y2 = gridCENTAUR%xyz(2,v2)
          z2 = gridCENTAUR%xyz(3,v2) 
          x3 = gridCENTAUR%xyz(1,v3)
          y3 = gridCENTAUR%xyz(2,v3)
          z3 = gridCENTAUR%xyz(3,v3)      
          
          xm = (x1 + x2 + x3)/3.0D0
          ym = (y1 + y2 + y3)/3.0D0
          zm = (z1 + z2 + z3)/3.0D0

          dx = xm - cofg(1,c2)
          dy = ym - cofg(2,c2)
          dz = zm - cofg(3,c2)

          cpx =   (y2 - y1)*(z3 - z1) - (z2 - z1)*(y3 - y1)
          cpy = -((x2 - x1)*(z3 - z1) - (z2 - z1)*(x3 - x1))
          cpz =   (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1) 

          IF ( cpx*dx + cpy*dy + cpz*dz < 0.0D0 ) THEN 
            gridCENTAUR%bQuad2v(1,l) = v1
            gridCENTAUR%bQuad2v(2,l) = v3
            gridCENTAUR%bQuad2v(3,l) = v2
          END IF ! cpx            
            
          v1 = gridCENTAUR%bQuad2v(1,l) ! get other triangular subface
          v2 = gridCENTAUR%bQuad2v(3,l)
          v3 = gridCENTAUR%bQuad2v(4,l)                   
          
          x1 = gridCENTAUR%xyz(1,v1)
          y1 = gridCENTAUR%xyz(2,v1)
          z1 = gridCENTAUR%xyz(3,v1)        
          x2 = gridCENTAUR%xyz(1,v2)
          y2 = gridCENTAUR%xyz(2,v2)
          z2 = gridCENTAUR%xyz(3,v2) 
          x3 = gridCENTAUR%xyz(1,v3)
          y3 = gridCENTAUR%xyz(2,v3)
          z3 = gridCENTAUR%xyz(3,v3)      
          
          xm = (x1 + x2 + x3)/3.0D0
          ym = (y1 + y2 + y3)/3.0D0
          zm = (z1 + z2 + z3)/3.0D0

          dx = xm - cofg(1,c2)
          dy = ym - cofg(2,c2)
          dz = zm - cofg(3,c2)

          cpx =   (y2 - y1)*(z3 - z1) - (z2 - z1)*(y3 - y1)
          cpy = -((x2 - x1)*(z3 - z1) - (z2 - z1)*(x3 - x1))
          cpz =   (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1) 

          IF ( cpx*dx + cpy*dy + cpz*dz < 0.0D0 ) THEN ! impossible?
            WRITE(*,*) 'ERROR - inconsistent boundary quad face.'
            STOP
          END IF ! cpx      
        END IF ! gridCOBALT
      END IF ! c1
    END DO ! j
  END DO ! i    

! =============================================================================
! Let user change patch types and add names
! =============================================================================

  WRITE(STDOUT,'(3X,A)') 'Enter patch types and names:'
  WRITE(STDOUT,'(5X,A)') 'CENTAUR boundary types:'
  WRITE(STDOUT,'(7X,A)') '100 - Inflow'
  WRITE(STDOUT,'(7X,A)') '200 - Outflow'
  WRITE(STDOUT,'(7X,A)') '300 - Solid wall'
  WRITE(STDOUT,'(7X,A)') '400 - Slip wall'
  WRITE(STDOUT,'(7X,A)') '500 - Symmtery plane'
  WRITE(STDOUT,'(7X,A)') '700 - Periodic boundary'
  
  DO i = 1,gridCENTAUR%nBounds
    WRITE(STDOUT,'(5X,A,1X,I3)') 'Boundary patch:',i
    WRITE(STDOUT,'(7X,A)') 'Enter patch type:' 
    READ(STDIN,*) gridCENTAUR%bInfo(1,i)
    WRITE(STDOUT,'(7X,A)') 'Enter patch name:' 
    READ(STDIN,*) gridCENTAUR%bName(i)    
  END DO ! i   

! ------------------------------------------------------------------------------
! Comment
! ------------------------------------------------------------------------------

  
  WRITE(STDOUT,'(1X,A)') 'Conversion completed.'
  
! ******************************************************************************
! End
! ******************************************************************************
  
END SUBROUTINE COBALT2CENTAUR
