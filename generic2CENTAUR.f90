! ******************************************************************************
!
! $Id: generic2CENTAUR.f90,v 1.3 2005/03/10 01:33:57 haselbac Exp $
!
! Filename: generic2CENTAUR.F90
!
! Purpose: Convert 2d grid into 3d grid in CENTAUR format through extrusion.
!
! Description: None.
!
! Input: None.
! 
! Output: None.
!
! Notes:
!   1. Boundaries are labelled as follows: First come the boundaries of the 
!      2d grid, and the last two boundaries are the additional ones due to 
!      the extrusion. 
!
! Author: Andreas Haselbacher
!
! Copyright: (c) 2001 by the University of Illinois
!
! RCS Revision history:
!
!   $Log: generic2CENTAUR.f90,v $
!   Revision 1.3  2005/03/10 01:33:57  haselbac
!   Added setting of case name
!
!   Revision 1.2  2005/03/09 03:57:07  haselbac
!   Changed formatted write statement
!
!   Revision 1.1  2005/03/05 17:25:53  haselbac
!   Initial revision
!
!   Revision 1.4  2004/10/28 17:00:29  haselbac
!   Added check for zero width
!
!   Revision 1.3  2004/07/19 18:55:33  haselbac
!   Added writing of info about grid
!
!   Revision 1.2  2004/07/17 22:10:59  haselbac
!   Changed coordinate specification
!
!   Revision 1.1  2003/03/07 15:26:34  haselbac
!   Initial revision
!
!
! ******************************************************************************

SUBROUTINE generic2CENTAUR

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
 
  INTEGER :: bTypeFrontBack,ib,iBQuads,iBTris,ie,iHexs,il,iPris,iq,it,iv, & 
             iv3d,nLayers,quadOffs,triOffs
  DOUBLE PRECISION :: zBack,zFront,zWidth
  CHARACTER*(MAX_STRING_LEN) :: bNameBack,bNameFront

! ******************************************************************************
! Start, get additional information for extrusion
! ******************************************************************************

  WRITE(STDOUT,'(/,1X,A)') 'Creating 3d-grid in CENTAUR format...'

  WRITE(STDOUT,'(3X,A)') 'Enter type of boundary for front/back boundaries:'
  WRITE(STDOUT,'(3X,A)') '  1 - Solid wall (no-slip): Type 300'
  WRITE(STDOUT,'(3X,A)') '  2 - Solid wall (slip): Type 400'
  WRITE(STDOUT,'(3X,A)') '  3 - Symmetry: Type 500'
  WRITE(STDOUT,'(3X,A)') '  4 - Periodic: Type 600'
  READ(STDIN,*) ib 

  SELECT CASE ( ib ) 
    CASE ( 1 ) 
      bTypeFrontBack = 300
      bNameFront     = 'Front no-slip wall'
      bNameBack      = 'Back no-slip wall'
    CASE ( 2 ) 
      bTypeFrontBack = 400
      bNameFront     = 'Front slip wall'
      bNameBack      = 'Back slip wall'
    CASE ( 3 ) 
      bTypeFrontBack = 500
      bNameFront     = 'Front symmetry'
      bNameBack      = 'Back symmetry' 
    CASE ( 4 ) 
      bTypeFrontBack = 600
      bNameFront     = 'Front periodic'
      bNameBack      = 'Back periodic'
    CASE DEFAULT
      CALL errorHandling(REACHED_DEFAULT,'generic2CENTAUR')
  END SELECT

  WRITE(STDOUT,'(3X,A)') 'Enter number of cell-layers to be created:'
  READ(STDIN,*) nLayers

  WRITE(STDOUT,'(3X,A)') 'Enter z-coordinate position of front boundary:'
  READ(STDIN,*) zFront

  WRITE(STDOUT,'(3X,A)') 'Enter width in z-coordinate direction:'
  READ(STDIN,*) zWidth 

  IF ( zWidth == 0.0D0 ) THEN 
    CALL errorHandling(ZEROWIDTH_ERROR)
  END IF ! zWdith

  zBack = zFront - ABS(zWidth)

! ******************************************************************************
! Set sizes and allocate memory for CENTAUR grid
! ******************************************************************************

  gridCENTAUR%nBounds = grid%nBounds + 2
  gridCENTAUR%nCells  = grid%nCells*nLayers
  gridCENTAUR%nHexs   = grid%nQuads*nLayers
  gridCENTAUR%nPris   = grid%nTris*nLayers
  gridCENTAUR%nTets   = 0
  gridCENTAUR%nPyrs   = 0
  gridCENTAUR%nVert   = grid%nVert*(nLayers+1)

  IF ( gridCENTAUR%nHexs > 0 ) THEN 
    ALLOCATE(gridCENTAUR%hex2v(8,gridCENTAUR%nHexs),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%hex2v',errorFlag)
    END IF ! errorFlag
  END IF ! gridCENTAUR

  IF ( gridCENTAUR%nPris > 0 ) THEN 
    ALLOCATE(gridCENTAUR%pri2v(6,gridCENTAUR%nPris),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%pri2v',errorFlag)
    END IF ! errorFlag
  END IF ! gridCENTAUR

  ALLOCATE(gridCENTAUR%xyz(3,gridCENTAUR%nVert),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%xyz',errorFlag)
  END IF ! errorFlag

  ALLOCATE(gridCENTAUR%bInfo(3,gridCENTAUR%nBounds),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%bInfo',errorFlag)
  END IF ! errorFlag

  ALLOCATE(gridCENTAUR%bName(gridCENTAUR%nBounds),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%bName',errorFlag)
  END IF ! errorFlag
  
  DO ib = 1,grid%nBounds
    gridCENTAUR%bName(ib)   = grid%bound(ib)%bName
    gridCENTAUR%bInfo(1,ib) = grid%bound(ib)%bType 
    gridCENTAUR%bInfo(2,ib) = 0 ! never any triangles on extruded boundaries

    IF ( ib == 1 ) THEN 
      gridCENTAUR%bInfo(3,ib) = grid%bound(ib)%nEdges*nLayers 
    ELSE 
      gridCENTAUR%bInfo(3,ib) = gridCENTAUR%bInfo(3,ib-1) & 
                              + grid%bound(ib)%nEdges*nLayers 
    END IF ! ib
  END DO ! ib
 
  gridCENTAUR%bName(grid%nBounds+1) = bNameBack 
  gridCENTAUR%bName(grid%nBounds+2) = bNameFront

  gridCENTAUR%bInfo(1,grid%nBounds+1) = bTypeFrontBack
  gridCENTAUR%bInfo(1,grid%nBounds+2) = bTypeFrontBack

  gridCENTAUR%bInfo(2,grid%nBounds+1) = gridCENTAUR%bInfo(2,grid%nBounds) & 
                                      + grid%nTris
  gridCENTAUR%bInfo(3,grid%nBounds+1) = gridCENTAUR%bInfo(3,grid%nBounds) & 
                                      + grid%nQuads
  gridCENTAUR%bInfo(2,grid%nBounds+2) = gridCENTAUR%bInfo(2,grid%nBounds+1) &
                                      + grid%nTris
  gridCENTAUR%bInfo(3,grid%nBounds+2) = gridCENTAUR%bInfo(3,grid%nBounds+1) &
                                      + grid%nQuads

  gridCENTAUR%nBTris  = gridCENTAUR%bInfo(2,gridCENTAUR%nBounds)
  gridCENTAUR%nBQuads = gridCENTAUR%bInfo(3,gridCENTAUR%nBounds)

  IF ( gridCENTAUR%nBTris > 0 ) THEN 
    ALLOCATE(gridCENTAUR%bTri2v(3,gridCENTAUR%nBTris),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%bTri2v',errorFlag)
    END IF ! errorFlag 
  END IF ! gridCENTAUR

  IF ( gridCENTAUR%nBQuads > 0 ) THEN 
    ALLOCATE(gridCENTAUR%bQuad2v(4,gridCENTAUR%nBQuads),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN 
      CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%bQuad2v',errorFlag)
    END IF ! errorFlag 
  END IF ! gridCENTAUR

! ******************************************************************************
! Create first cell layer 
! ******************************************************************************

  iPris = 0
  iHexs = 0

  DO it = 1,grid%nTris
    iPris = iPris + 1
    
    gridCENTAUR%pri2v(1:3,iPris) = grid%tri2v(1:3,it)
    gridCENTAUR%pri2v(4:6,iPris) = gridCENTAUR%pri2v(1:3,iPris) + grid%nVert
  END DO ! it

  DO iq = 1,grid%nQuads
    iHexs = iHexs + 1
    
    gridCENTAUR%hex2v(1:4,iHexs) = grid%quad2v(1:4,iq)
    gridCENTAUR%hex2v(5:8,iHexs) = gridCENTAUR%hex2v(1:4,iHexs) + grid%nVert    
  END DO ! iq

! ******************************************************************************
! Create second to final cell layer 
! ******************************************************************************

  DO il = 2,nLayers
    DO it = 1,grid%nTris
      iPris = iPris + 1

      gridCENTAUR%pri2v(1:3,iPris) = gridCENTAUR%pri2v(4:6,iPris-grid%nTris) 
      gridCENTAUR%pri2v(4:6,iPris) = gridCENTAUR%pri2v(1:3,iPris) + grid%nVert
    END DO ! it

    DO iq = 1,grid%nQuads
      iHexs = iHexs + 1

      gridCENTAUR%hex2v(1:4,iHexs) = gridCENTAUR%hex2v(5:8,iHexs-grid%nQuads) 
      gridCENTAUR%hex2v(5:8,iHexs) = gridCENTAUR%hex2v(1:4,iHexs) + grid%nVert
    END DO ! iq 
  END DO ! il

! ******************************************************************************
! Create coordinates
! ******************************************************************************

  DO il = 1,nLayers+1
    DO iv = 1,grid%nVert
      iv3d = iv + (il-1)*grid%nVert 

      gridCENTAUR%xyz(1,iv3d) = grid%xy(1,iv)
      gridCENTAUR%xyz(2,iv3d) = grid%xy(2,iv)
      gridCENTAUR%xyz(3,iv3d) = zBack + (il-1)/REAL(nLayers)*(zFront-zBack) 
    END DO ! iv
  END DO ! il 

! ******************************************************************************
! Create boundary faces for extruded boundaries 
! ******************************************************************************

  iBTris  = 0
  iBQuads = 0

  DO ib = 1,grid%nBounds
  
! =============================================================================
!   First layer of faces 
! ============================================================================= 
    
    DO ie = 1,grid%bound(ib)%nEdges
      iBQuads = iBQuads + 1
      
      gridCENTAUR%bQuad2v(1,iBQuads) = grid%bound(ib)%e2v(1,ie)
      gridCENTAUR%bQuad2v(2,iBQuads) = grid%bound(ib)%e2v(2,ie) 
      gridCENTAUR%bQuad2v(3,iBQuads) = grid%bound(ib)%e2v(2,ie) + grid%nVert
      gridCENTAUR%bQuad2v(4,iBQuads) = grid%bound(ib)%e2v(1,ie) + grid%nVert    
    END DO ! ie

! =============================================================================
!   Second to last layer of faces 
! ============================================================================= 

    DO il = 2,nLayers
      DO ie = 1,grid%bound(ib)%nEdges
        iBQuads = iBQuads + 1

        gridCENTAUR%bQuad2v(1,iBQuads) = & 
          gridCENTAUR%bQuad2v(1,iBQuads-grid%bound(ib)%nEdges) + grid%nVert
        gridCENTAUR%bQuad2v(2,iBQuads) = & 
          gridCENTAUR%bQuad2v(2,iBQuads-grid%bound(ib)%nEdges) + grid%nVert
        gridCENTAUR%bQuad2v(3,iBQuads) = & 
          gridCENTAUR%bQuad2v(3,iBQuads-grid%bound(ib)%nEdges) + grid%nVert
        gridCENTAUR%bQuad2v(4,iBQuads) = & 
          gridCENTAUR%bQuad2v(4,iBQuads-grid%bound(ib)%nEdges) + grid%nVert
      END DO ! ie
    END DO ! il
  END DO ! ib

! ******************************************************************************
! Create boundary faces for front and back boundaries 
! ******************************************************************************

! ==============================================================================
! Front face - NOTE change of orientation
! ==============================================================================

  IF ( grid%nTris > 0 ) THEN 
    DO it = 1,grid%nTris
      iBTris = iBTris + 1

      gridCENTAUR%bTri2v(1,iBTris) = grid%tri2v(1,it)
      gridCENTAUR%bTri2v(2,iBTris) = grid%tri2v(3,it)
      gridCENTAUR%bTri2v(3,iBTris) = grid%tri2v(2,it)
    END DO ! it
  END IF ! grid

  IF ( grid%nQuads > 0 ) THEN 
    DO iq = 1,grid%nQuads
      iBQuads = iBQuads + 1

      gridCENTAUR%bQuad2v(1,iBQuads) = grid%quad2v(1,iq)
      gridCENTAUR%bQuad2v(2,iBQuads) = grid%quad2v(4,iq)
      gridCENTAUR%bQuad2v(3,iBQuads) = grid%quad2v(3,iq)
      gridCENTAUR%bQuad2v(4,iBQuads) = grid%quad2v(2,iq)
    END DO ! iq
  END IF ! grid
  
! ==============================================================================
! Back face - NOTE change of orientation
! ==============================================================================

  IF ( grid%nTris > 0 ) THEN 
    DO it = 1,grid%nTris
      iBTris = iBTris + 1

      gridCENTAUR%bTri2v(1,iBTris) = grid%tri2v(1,it) + nLayers*grid%nVert
      gridCENTAUR%bTri2v(2,iBTris) = grid%tri2v(2,it) + nLayers*grid%nVert
      gridCENTAUR%bTri2v(3,iBTris) = grid%tri2v(3,it) + nLayers*grid%nVert
    END DO ! it
  END IF ! grid

  IF ( grid%nQuads > 0 ) THEN 
    DO iq = 1,grid%nQuads
      iBQuads = iBQuads + 1

      gridCENTAUR%bQuad2v(1,iBQuads) = grid%quad2v(1,iq) + nLayers*grid%nVert
      gridCENTAUR%bQuad2v(2,iBQuads) = grid%quad2v(2,iq) + nLayers*grid%nVert
      gridCENTAUR%bQuad2v(3,iBQuads) = grid%quad2v(3,iq) + nLayers*grid%nVert
      gridCENTAUR%bQuad2v(4,iBQuads) = grid%quad2v(4,iq) + nLayers*grid%nVert
    END DO ! iq
  END IF ! grid      
    
  WRITE(STDOUT,'(1X,A)') 'Grid created successfully.' 
  
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
  DO ib = 1,gridCENTAUR%nBounds
    IF ( ib > 1 ) THEN
      triOffs  = gridCENTAUR%bInfo(2,ib-1)
      quadOffs = gridCENTAUR%bInfo(3,ib-1)           
    ELSE
      triOffs  = 0
      quadOffs = 0      
    END IF ! ib    
  
    WRITE(STDOUT,'(3X,I2,1X,A20,1X,I3,2(1X,I7)))') & 
      ib,TRIM(gridCENTAUR%bName(ib)), &
      gridCENTAUR%bInfo(1,ib), & 
      gridCENTAUR%bInfo(2,ib) - triOffs, & 
      gridCENTAUR%bInfo(3,ib) - quadOffs                      
  END DO ! ib

! ******************************************************************************
! Set case name
! ******************************************************************************
 
  gridCENTAUR%title = 'conflu' 
     
! ******************************************************************************
! End
! ******************************************************************************
  
END SUBROUTINE generic2CENTAUR
