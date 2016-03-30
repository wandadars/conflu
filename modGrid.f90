! ******************************************************************************
!
! $Id: modGrid.f90,v 1.3 2005/03/12 03:14:49 haselbac Exp $
!
! Filename: modGrid.F90
!
! Purpose: Collect type definitions for grid data structures
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
!   $Log: modGrid.f90,v $
!   Revision 1.3  2005/03/12 03:14:49  haselbac
!   Filled GAMBIT grid type
!
!   Revision 1.2  2005/03/10 22:06:30  haselbac
!   Added GAMBIT type
!
!   Revision 1.1  2005/03/05 17:25:53  haselbac
!   Initial revision
!
!   Revision 1.3  2004/12/27 15:32:00  haselbac
!   Added PLOT3D type and vMatch array
!
!   Revision 1.2  2004/07/17 22:12:09  haselbac
!   Added 3d boundary and volume grid types for CENTAUR conversion
!
!   Revision 1.1  2003/03/07 15:26:34  haselbac
!   Initial revision
!
! ******************************************************************************

MODULE modGrid

  USE modGlobals

  IMPLICIT NONE
  
  SAVE
  
! ******************************************************************************
! Generic to grids from all grid generators
! ******************************************************************************
  
! ==============================================================================
! Boundary type 2d
! ==============================================================================  
  
  TYPE bound_t
    INTEGER :: bType,nEdges,nVert
    INTEGER, DIMENSION(:), POINTER :: e2c,v
    INTEGER, DIMENSION(:,:), POINTER :: e2v
    CHARACTER*(MAX_STRING_LEN) :: bName
  END TYPE bound_t

! ==============================================================================
! Boundary type 3d
! ==============================================================================  
  
  TYPE bound3d_t
    LOGICAL :: selectFlag
    INTEGER :: nBEdges,nEdges,nQuads,nTris,nVert   
    INTEGER, DIMENSION(:), POINTER :: be2c,v,vMatch 
    INTEGER, DIMENSION(:,:), POINTER :: be2v,e2c,e2v,quad2v,tri2v
  END TYPE bound3d_t

! ==============================================================================
! Grid type 2d
! ==============================================================================  
 
  TYPE grid_t
    INTEGER :: nBounds,nCells,nQuads,nTris,nVert
    INTEGER, DIMENSION(:,:), POINTER :: quad2v,tri2v
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: xy
    TYPE(bound_t), DIMENSION(:), POINTER :: bound
  END TYPE grid_t

  TYPE(grid_t) :: grid

! ==============================================================================
! Grid type 3d
! ==============================================================================  
 
  TYPE grid3d_t
    INTEGER :: nBounds
    INTEGER, DIMENSION(:), POINTER :: boundMap    
    TYPE(bound3d_t), DIMENSION(:), POINTER :: bound
  END TYPE grid3d_t

  TYPE(grid3d_t) :: grid3d

! ******************************************************************************
! Specific to grids from COBALT solver (Matthew Grismer, WPAFB) 
! ******************************************************************************

  TYPE gridCOBALT_t
    INTEGER :: nDim,nFaces,nFacesPerCellMax,nQuads,nVertPerFaceMax,nTris, & 
               nZones,nPatches,nMappings
    INTEGER, DIMENSION(:), POINTER :: nvpf
    INTEGER, DIMENSION(:,:), POINTER :: f2c,f2v,patch2bc
  END TYPE gridCOBALT_t

  TYPE(gridCOBALT_t) :: gridCOBALT

! ******************************************************************************
! Specific to grids from CENTAUR (Yannis Kallinderis, CentaurSoft)
! ******************************************************************************

  TYPE gridCENTAUR_t 
    INTEGER :: nBounds,nBQuads,nBTris,nCells,nHexs,nPris,nPyrs,nTets,nVert
    INTEGER, DIMENSION(:,:), POINTER :: bInfo,bTri2v,bQuad2v,hex2v,pri2v, & 
                                        pyr2v,tet2v 
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: xyz
    CHARACTER*(MAX_STRING_LEN) :: title
    CHARACTER*(MAX_STRING_LEN), DIMENSION(:), POINTER :: bName    
  END TYPE gridCENTAUR_t

  TYPE(gridCENTAUR_t) :: gridCENTAUR

! ******************************************************************************
! Specific to grids from GAMBIT (Fluent) 
! ******************************************************************************

  TYPE gridGAMBIT_t
    CHARACTER*(MAX_STRING_LEN) :: title   
    INTEGER :: NELGP,NDP,NFLAGS,NTYPE
    INTEGER :: nCells,nPatches,nQuads,nTris,nVert
    INTEGER, DIMENSION(:), POINTER :: ct
    INTEGER, DIMENSION(:,:), POINTER :: c2v
    TYPE(patchGAMBIT_t), DIMENSION(:), POINTER :: patches
  END TYPE gridGAMBIT_t

  TYPE patchGAMBIT_t
    CHARACTER*(MAX_STRING_LEN) :: patchName   
    INTEGER :: bcType,nFaces
    INTEGER, DIMENSION(:), POINTER :: bf2ct,bf2cgi,bf2fli
  END TYPE patchGAMBIT_t

  TYPE(gridGAMBIT_t) :: gridGAMBIT

! ******************************************************************************
! Specific to PLOT3D grids
! ******************************************************************************
                                                                                                                                             
  TYPE gridPLOT3D_t
    INTEGER :: nBounds,nCells,nVert,nx,ny,nz
    DOUBLE PRECISION, DIMENSION(:,:,:,:), POINTER :: xyz
  END TYPE gridPLOT3D_t
                                                                                                                                             
  TYPE(gridPLOT3D_t) :: gridPLOT3D

! ******************************************************************************
! Mapping from faces to vertices for each cell type: required for conversion
! of COBALT to ROCFLU format. NOTE the vertices are ordered anticlockwise as 
! to be consistent with outward normal vectors  
! ******************************************************************************  
 
  INTEGER, DIMENSION(4,4) :: f2vTet = & 
    RESHAPE((/1,2,3,0,2,4,3,0,1,3,4,0,1,4,2,0/), (/4,4/))
  INTEGER, DIMENSION(4,6) :: f2vHex = & 
    RESHAPE((/1,4,3,2,1,2,6,5,2,3,7,6,3,4,8,7,1,5,8,4,5,6,7,8/), (/4,6/))
  INTEGER, DIMENSION(4,5) :: f2vPri = & 
    RESHAPE((/1,3,2,0,1,2,5,4,2,3,6,5,1,4,6,3,4,5,6,0/), (/4,5/))
  INTEGER, DIMENSION(4,5) :: f2vPyr = & 
    RESHAPE((/1,4,3,2,1,2,5,0,2,3,5,0,3,4,5,0,1,5,4,0/), (/4,5/))    
        
! ******************************************************************************
! End
! ******************************************************************************

END MODULE modGrid

