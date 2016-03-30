! ******************************************************************************
!
! $Id: readGridGAMBIT2d.f90,v 1.2 2005/03/12 03:18:08 haselbac Exp $
!
! Filename: readGridGAMBIT2d.F90
!
! Purpose: Read 2d grid in GAMBIT neutral file format.
!
! Description: None.
!
! Input: None.
! 
! Output: None.
!
! Notes: 
!   1. Restricted to two-dimensional GAMBIT grids.
!
! Author: Andreas Haselbacher and Fady Najjar
!
! Copyright: (c) 2005 by the University of Illinois
!
! RCS Revision history:
!
!   $Log: readGridGAMBIT2d.f90,v $
!   Revision 1.2  2005/03/12 03:18:08  haselbac
!   Filled routine, minimal testing so far
!
!   Revision 1.1  2005/03/10 22:06:10  haselbac
!   Initial revision
!
! ******************************************************************************

SUBROUTINE readGridGAMBIT2d

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
  
  CHARACTER*(MAX_STRING_LEN) :: dummyString,iFileName,sectionString,&
                                versionString  
  INTEGER :: dummyInteger,ib,icg,icl,ict,icv,iFile,ifl,ifl2,iPatch,ITYPE, &
             ivg,ivl,ivl1,ivl2,loopCounter,NENTRY
  TYPE(patchGAMBIT_t), POINTER :: pPatchGAMBIT

! ==============================================================================
! Element types
! ==============================================================================

  INTEGER, PARAMETER :: GAMBIT_NTYPE_EDGE = 1, & 
                        GAMBIT_NTYPE_QUAD = 2, & 
                        GAMBIT_NTYPE_TRI  = 3, & 
                        GAMBIT_NTYPE_HEX  = 4, & 
                        GAMBIT_NTYPE_PRI  = 5, &
                        GAMBIT_NTYPE_TET  = 6, & 
                        GAMBIT_NTYPE_PYR  = 7
 
! ==============================================================================
! Mapping of faces to vertices for GAMBIT
! ==============================================================================
                                                       
  INTEGER, DIMENSION(2,3), PARAMETER :: f2vTriGAMBIT = &
    RESHAPE((/1,2,2,3,3,1/), (/2,3/))
  INTEGER, DIMENSION(2,4), PARAMETER :: f2vQuadGAMBIT = &
    RESHAPE((/1,2,2,3,3,4,4,1/), (/2,4/))
      
! ******************************************************************************
! Start, open file
! ******************************************************************************
  
  WRITE(STDOUT,'(1X,A)') 'Reading GAMBIT2D neutral grid file....'
  WRITE(STDOUT,'(3X,A)') 'Enter file name:'
  READ(STDIN,*) iFileName

  iFile = FILE_UNIT_GRID_INPUT

  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)   
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(FILE_OPEN_ERROR,iFileName,errorFlag)
  END IF ! errorFlag 

! ******************************************************************************
! Initialize variables
! ******************************************************************************

  gridGAMBIT%NDP   = 0
  gridGAMBIT%NTYPE = 0

  gridGAMBIT%nQuads   = 0
  gridGAMBIT%nTris    = 0
  gridGAMBIT%nPatches = 0
  gridGAMBIT%nVert    = 0

  iPatch = 0

! ******************************************************************************
! Read header
! ******************************************************************************

  READ(iFile,'(2(A20))') dummyString,versionString
  READ(iFile,*) dummyString    
  READ(iFile,*) gridGAMBIT%title
  READ(iFile,*) dummyString 
  READ(iFile,*) dummyString
  READ(iFile,*) dummyString
  READ(iFile,*) gridGAMBIT%nVert,gridGAMBIT%nCells,dummyInteger, &
                gridGAMBIT%nPatches,dummyInteger,dummyInteger
  READ(iFile,*) dummyString

  IF ( TRIM(dummyString) /= "ENDOFSECTION" ) THEN 
    CALL errorHandling(STRING_INVALID)
  END IF ! TRIM(sectionString)     

! ******************************************************************************
! Allocate memory
! ******************************************************************************

  grid%nVert  = gridGAMBIT%nVert
  grid%nCells = gridGAMBIT%nCells  

  ALLOCATE(grid%xy(2,grid%nVert),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'grid%xy',errorFlag)
  END IF ! errorFlag

  ALLOCATE(gridGAMBIT%ct(gridGAMBIT%nCells),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'gridGAMBIT%ct',errorFlag)
  END IF ! errorFlag

  ALLOCATE(gridGAMBIT%c2v(4,gridGAMBIT%nCells),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'gridGAMBIT%c2v',errorFlag)
  END IF ! errorFlag

  ALLOCATE(gridGAMBIT%patches(gridGAMBIT%nPatches),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(ALLOCATE_ERROR,'gridGAMBIT%patches',errorFlag)
  END IF ! errorFlag

! ******************************************************************************
! Read sections
! ******************************************************************************
  
  loopCounter = 0

  DO ! set up infinite loop
    loopCounter = loopCounter + 1

! ==============================================================================
!   Read section string and take appropriate action
! ==============================================================================
  
    READ(iFile,'(A32)',IOSTAT=errorFlag,END=1) sectionString

! ==============================================================================
!   Coordinates
! ==============================================================================

    IF ( ADJUSTL(TRIM(sectionString)) == & 
         "NODAL COORDINATES"//" "//ADJUSTL(TRIM(versionString)) ) THEN           
      WRITE(STDOUT,'(3X,A)') 'Reading GAMBIT coordinate section...'         

      DO ivl = 1,gridGAMBIT%nVert
        READ(iFile,*) ivg
        BACKSPACE(iFile)

        READ(iFile,*) dummyInteger,grid%xy(1:2,ivg)
      END DO ! ivl

! ==============================================================================
!   Element connectivity
! ==============================================================================
      
    ELSE IF ( ADJUSTL(TRIM(sectionString)) == & 
              "ELEMENTS/CELLS"//" "//ADJUSTL(TRIM(versionString)) ) THEN 
      WRITE(STDOUT,'(3X,A)') 'Reading GAMBIT element connectivity...'

      DO icl = 1,gridGAMBIT%nCells
        READ(iFile,*) icg
        BACKSPACE(iFile)

        READ(iFile,*) dummyString,gridGAMBIT%NTYPE,gridGAMBIT%NDP
        BACKSPACE(iFile)

        gridGAMBIT%ct(icg) = gridGAMBIT%NTYPE

! ------------------------------------------------------------------------------
!       Read element connectivity
! ------------------------------------------------------------------------------
          
        SELECT CASE ( gridGAMBIT%NTYPE ) 
            
! ------- Edge -----------------------------------------------------------------   
            
          CASE ( GAMBIT_NTYPE_EDGE ) 
            CALL errorHandling(NTYPE_INVALID_ERROR)             

! ------- Quadrilateral --------------------------------------------------------        

          CASE ( GAMBIT_NTYPE_QUAD )            
            IF ( gridGAMBIT%NDP == 4 ) THEN
              gridGAMBIT%nQuads = gridGAMBIT%nQuads + 1

              READ(iFile,*) dummyString,dummyString,dummyString, & 
                            gridGAMBIT%c2v(1:4,icg)                                                     
            ELSE 
              CALL errorHandling(NDP_INVALID_ERROR)
            END IF ! gridGAMBIT%NDP

! ------- Triangle -------------------------------------------------------------      

          CASE ( GAMBIT_NTYPE_TRI )
            CALL errorHandling(NTYPE_INVALID_ERROR)

! ------- Hexahedron -----------------------------------------------------------           

          CASE ( GAMBIT_NTYPE_HEX ) 
             CALL errorHandling(NTYPE_INVALID_ERROR)            

! ------- Prism ----------------------------------------------------------------            

          CASE ( GAMBIT_NTYPE_PRI ) 
            CALL errorHandling(NTYPE_INVALID_ERROR)                        

! ------- Tetrahedron ----------------------------------------------------------             

          CASE ( GAMBIT_NTYPE_TET ) 
            CALL errorHandling(NTYPE_INVALID_ERROR)                                

! ------- Pyramid --------------------------------------------------------------             

          CASE ( GAMBIT_NTYPE_PYR ) 
            CALL errorHandling(NTYPE_INVALID_ERROR)                  
             
! ------- Unknown --------------------------------------------------------------              
                            
          CASE DEFAULT
            CALL errorHandling(REACHED_DEFAULT)   
        END SELECT ! gridGAMBIT%NTYPE                                 
      END DO ! icl

! ==============================================================================
!   Element group
! ==============================================================================

    ELSE IF ( ADJUSTL(TRIM(sectionString)) == & 
              "ELEMENT GROUP"//" "//ADJUSTL(TRIM(versionString)) ) THEN        
      WRITE(STDOUT,'(3X,A)') 'Reading GAMBIT element group information...'

      READ(iFile,*) dummyString,dummyInteger,dummyString,gridGAMBIT%NELGP, & 
                    dummyString,dummyInteger,dummyString,gridGAMBIT%NFLAGS 

      READ(iFile,*) dummyString          

      READ(iFile,'(10(I8))') (dummyInteger,ifl=1,gridGAMBIT%NFLAGS)
      READ(iFile,'(10(I8))') (dummyInteger,icl=1,gridGAMBIT%NELGP)

! ==============================================================================
!   Boundary information
! ==============================================================================

    ELSE IF ( ADJUSTL(TRIM(sectionString)) == & 
              "BOUNDARY CONDITIONS"//" "//ADJUSTL(TRIM(versionString)) ) THEN 
      WRITE(STDOUT,'(3X,A)') 'Reading GAMBIT boundary condition information...'

      iPatch = iPatch + 1

      pPatchGAMBIT => gridGAMBIT%patches(iPatch)   

      READ(iFile,'(A32,2(I10))') pPatchGAMBIT%patchName,ITYPE,NENTRY

! ------------------------------------------------------------------------------
!     Read face or vertex data. NOTE in either case, values are not read. 
! ------------------------------------------------------------------------------

      SELECT CASE ( ITYPE )

! ----- Face data --------------------------------------------------------------  

        CASE ( 1 ) ! Face data              
          pPatchGAMBIT%nFaces = NENTRY

          ALLOCATE(pPatchGAMBIT%bf2cgi(pPatchGAMBIT%nFaces),STAT=errorFlag)
          IF ( errorFlag /= NO_ERROR ) THEN 
            CALL errorHandling(ALLOCATE_ERROR,'pPatchGAMBIT%bf2cgi',errorFlag)
          END IF ! errorFlag

          ALLOCATE(pPatchGAMBIT%bf2ct(pPatchGAMBIT%nFaces),STAT=errorFlag)
          IF ( errorFlag /= NO_ERROR ) THEN 
            CALL errorHandling(ALLOCATE_ERROR,'pPatchGAMBIT%bf2ct',errorFlag)
          END IF ! errorFlag

          ALLOCATE(pPatchGAMBIT%bf2fli(pPatchGAMBIT%nFaces),STAT=errorFlag)
          IF ( errorFlag /= NO_ERROR ) THEN 
            CALL errorHandling(ALLOCATE_ERROR,'pPatchGAMBIT%bf2fli',errorFlag)
          END IF ! errorFlag             

          DO ifl = 1,pPatchGAMBIT%nFaces
            READ(iFile,*) pPatchGAMBIT%bf2cgi(ifl), & 
                          pPatchGAMBIT%bf2ct(ifl), & 
                          pPatchGAMBIT%bf2fli(ifl)                              
          END DO ! ifl

! ----- Vertex data ------------------------------------------------------------

        CASE ( 0 ) ! Vertex data
          DO ivl = 1,NENTRY
            READ(iFile,*) dummyInteger
          END DO ! ivl

! ----- Unknown ----------------------------------------------------------------             

        CASE DEFAULT 
          CALL errorHandling(REACHED_DEFAULT)          
      END SELECT ! ITYPE 
          
! ==============================================================================
!   Unknown section string
! ==============================================================================
      
    ELSE
      CALL errorHandling(REACHED_DEFAULT)      
    END IF ! TRIM(sectionString)

! ==============================================================================
!   Read end of section string 
! ==============================================================================

    READ(iFile,*) dummyString
      
    IF ( TRIM(dummyString) /= "ENDOFSECTION" ) THEN 
      CALL errorHandling(STRING_INVALID)
    END IF ! TRIM(sectionString)

! ==============================================================================
!   Guard against infinite loop - might be unnecessary because of read errors?
! ==============================================================================

    IF ( loopCounter >= 1E6 ) THEN
      CALL errorHandling(INFINITE_LOOP)
    END IF ! loopCounter
  END DO ! <empty>

! ******************************************************************************
! EOF condition. NOTE assume that this is ok, because GAMBIT grid files do not
! have a dedicated EOF marker, and hence need to use READ statement to detect
! EOF condition.
! ******************************************************************************

1 WRITE(STDOUT,*) '*** WARNING *** Encountered EOF.'

  WRITE(STDOUT,'(1X,A)') 'Reading GAMBIT neutral grid file done.'

! ******************************************************************************
! Print grid statistics
! ******************************************************************************

  WRITE(STDOUT,'(/,1X,A)')     'GAMBIT Grid Statistics:'
  WRITE(STDOUT,'(3X,A,2X,I9)') 'Vertices:       ',gridGAMBIT%nVert
  WRITE(STDOUT,'(3X,A,2X,I9)') 'Cells:          ',gridGAMBIT%nCells
  WRITE(STDOUT,'(5X,A,I9)')    'Triangles:      ',gridGAMBIT%nTris  
  WRITE(STDOUT,'(5X,A,I9)')    'Quadrilaterals: ',gridGAMBIT%nQuads

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(iFile,IOSTAT=errorFlag)  
  IF ( errorFlag /= NO_ERROR ) THEN 
    CALL errorHandling(FILE_CLOSE_ERROR,iFileName,errorFlag)
  END IF ! global%error  

! ******************************************************************************
! Convert to generic format
! ******************************************************************************

! =============================================================================
! Connectivity
! =============================================================================

  grid%nTris  = gridGAMBIT%nTris
  grid%nQuads = gridGAMBIT%nQuads

  ALLOCATE(grid%tri2v(3,grid%nTris),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'grid%tri2v',errorFlag)
  END IF ! errorFlag  

  ALLOCATE(grid%quad2v(4,grid%nQuads),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'grid%quad2v',errorFlag)
  END IF ! errorFlag  

  grid%nTris  = 0
  grid%nQuads = 0

  DO icg = 1,gridGAMBIT%nCells
    SELECT CASE ( gridGAMBIT%ct(icg) )
      CASE ( GAMBIT_NTYPE_TRI ) 
        grid%nTris = grid%nTris + 1
      
        grid%tri2v(1,grid%nTris) = gridGAMBIT%c2v(1,icg)
        grid%tri2v(2,grid%nTris) = gridGAMBIT%c2v(2,icg)
        grid%tri2v(3,grid%nTris) = gridGAMBIT%c2v(3,icg)                
      CASE ( GAMBIT_NTYPE_QUAD )
        grid%nQuads = grid%nQuads + 1      
      
        grid%quad2v(1,grid%nQuads) = gridGAMBIT%c2v(1,icg)
        grid%quad2v(2,grid%nQuads) = gridGAMBIT%c2v(2,icg)
        grid%quad2v(3,grid%nQuads) = gridGAMBIT%c2v(3,icg)
        grid%quad2v(4,grid%nQuads) = gridGAMBIT%c2v(4,icg)            
      CASE DEFAULT
        CALL errorHandling(REACHED_DEFAULT)
    END SELECT ! gridGAMBIT%ct 
  END DO ! icg

! =============================================================================
! Boundary
! =============================================================================

  grid%nBounds = gridGAMBIT%nPatches
  
  ALLOCATE(grid%bound(grid%nBounds),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'grid%bound',errorFlag)
  END IF ! errorFlag  

  DO ib = 1,grid%nBounds
    pPatchGAMBIT => gridGAMBIT%patches(ib)    
  
    grid%bound(ib)%nEdges = pPatchGAMBIT%nFaces
    grid%bound(ib)%bName  = ADJUSTL(TRIM(pPatchGAMBIT%patchName))
        
    ALLOCATE(grid%bound(ib)%e2v(2,grid%bound(ib)%nEdges),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN
      CALL errorHandling(ALLOCATE_ERROR,'grid%bound%e2v',errorFlag)
    END IF ! errorFlag 
  END DO ! ib  

  DO ib = 1,grid%nBounds
    pPatchGAMBIT => gridGAMBIT%patches(ib)    
  
    DO ifl = 1,grid%bound(ib)%nEdges
      ict  = pPatchGAMBIT%bf2ct(ifl)
      icg  = pPatchGAMBIT%bf2cgi(ifl)
      ifl2 = pPatchGAMBIT%bf2fli(ifl)
      
      SELECT CASE ( ict ) 
        CASE ( GAMBIT_NTYPE_QUAD )               
          ivl1 = f2vQuadGAMBIT(1,ifl2)
          ivl2 = f2vQuadGAMBIT(2,ifl2)
                  
          grid%bound(ib)%e2v(1,ifl) = gridGAMBIT%c2v(ivl1,icg)
          grid%bound(ib)%e2v(2,ifl) = gridGAMBIT%c2v(ivl2,icg)            
        CASE DEFAULT
          CALL errorHandling(REACHED_DEFAULT) 
      END SELECT ! ict
    END DO ! ifl
  END DO ! ib   

! ==============================================================================
! Assign boundary conditions to existing boundaries 
! ==============================================================================

  WRITE(STDOUT,'(/,1X,A)') 'Assign boundary numbers from following selection:'
  WRITE(STDOUT,'(3X,A)') '100 - Inflow'
  WRITE(STDOUT,'(3X,A)') '200 - Outflow'
  WRITE(STDOUT,'(3X,A)') '300 - No-slip wall'
  WRITE(STDOUT,'(3X,A)') '400 - Slip wall'
  WRITE(STDOUT,'(3X,A)') '500 - Farfield'

  DO ib = 1,grid%nBounds
    WRITE(STDOUT,'(3X,A,1X,I2)') 'Boundary:',ib  
    WRITE(STDOUT,'(5X,A,1X,A)') 'GAMBIT boundary name:',grid%bound(ib)%bName

    WRITE(STDOUT,'(5X,A)') 'Enter boundary number:'
    READ(STDIN,*) grid%bound(ib)%bType
  END DO ! ib
  
! ******************************************************************************
! End
! ******************************************************************************
  
END SUBROUTINE readGridGAMBIT2d
