! ******************************************************************************
!
! $Id: PLOT3D2CENTAUR.f90,v 1.1 2005/03/05 17:25:53 haselbac Exp $
!
! Filename: PLOT3D2CENTAUR.F90
!
! Purpose: Convert PLOT3D grid format into CENTAUR grid format.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes: 
!  1. So far tested only for one geometry, so may still have bugs, especially 
!     in removing wake cut.
!  2. Not general, works only for grids with one boundary condition for each 
!     block face. 
!
! Author: Andreas Haselbacher
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! RCS Revision history:
!
! $Log: PLOT3D2CENTAUR.f90,v $
! Revision 1.1  2005/03/05 17:25:53  haselbac
! Initial revision
!
! Revision 1.2  2005/02/18 20:37:26  haselbac
! Added comments, deleted superfluous RCS log entry
!
! Revision 1.1  2004/12/27 15:34:48  haselbac
! Initial revision
!
! ******************************************************************************

SUBROUTINE PLOT3D2CENTAUR

  USE modError
  USE modGlobals
  USE modGrid
  USE modHashTable
  USE modSortSearch

  IMPLICIT NONE

  INTERFACE
    INTEGER FUNCTION ijk2l(i,j,k,ni,nj)
      INTEGER, INTENT(IN) :: i,j,k,ni,nj
    END FUNCTION
  END INTERFACE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Parameters
! ==============================================================================

! ==============================================================================
! Local variables
! ==============================================================================

  INTEGER :: cntr,faceType,ib,ic,iq,iq2,iv,iv1,iv2,ix,iy,iz,key,l,nFaces, & 
             nMatches,nVertMatches,nVertNew,nx,ny,nz,v1,v2
  INTEGER :: fv(4)
  INTEGER, DIMENSION(:), ALLOCATABLE :: vMatch,vMatchInv
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: f2v,vmapInv 
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: vmap
  DOUBLE PRECISION :: dx,dy,dz,toler

! ******************************************************************************
! Start
! ******************************************************************************

  WRITE(STDOUT,'(/,1X,A)') 'Converting from PLOT3D to CENTAUR format...'


  gridCENTAUR%nVert = gridPLOT3d%nVert

! ******************************************************************************
! Determine vertex map
! ******************************************************************************

  WRITE(STDOUT,'(3X,A)') 'Determining vertex mapping...'

  nx = gridPLOT3D%nx 
  ny = gridPLOT3D%ny
  nz = gridPLOT3D%nz

  ALLOCATE(vmap(nx,ny,nz),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'vmap',errorFlag)
  END IF ! errorFlag
  
  ALLOCATE(vmapInv(3,gridCENTAUR%nVert),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'vmapInv',errorFlag)
  END IF ! errorFlag  
  
  DO iz = 1,nz
    DO iy = 1,ny
      DO ix = 1,nx      
        l = ijk2l(ix,iy,iz,nx,ny)
        
        vmap(ix,iy,iz) = l
         
        vmapInv(1,l) = ix
        vmapInv(2,l) = iy
        vmapInv(3,l) = iz        
      END DO ! ix
    END DO ! iy
  END DO ! iz  

! ******************************************************************************
! Copy vertices
! ******************************************************************************
  
  ALLOCATE(gridCENTAUR%xyz(3,gridCENTAUR%nVert),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%xyz',errorFlag)
  END IF ! errorFlag  

  DO iz = 1,nz
    DO iy = 1,ny
      DO ix = 1,nx      
        iv = vmap(ix,iy,iz) 
        
        gridCENTAUR%xyz(1,iv) = gridPLOT3D%xyz(1,ix,iy,iz)
        gridCENTAUR%xyz(2,iv) = gridPLOT3D%xyz(2,ix,iy,iz)
        gridCENTAUR%xyz(3,iv) = gridPLOT3D%xyz(3,ix,iy,iz)                
      END DO ! ix
    END DO ! iy
  END DO ! iz  

! ******************************************************************************
! Determine cell connectivity
! ******************************************************************************

  WRITE(STDOUT,'(3X,A)') 'Determining cell mapping...'
  
  gridCENTAUR%nTets = 0
  gridCENTAUR%nHexs = gridPLOT3D%nCells
  gridCENTAUR%nPris = 0
  gridCENTAUR%nPyrs = 0

  ALLOCATE(gridCENTAUR%hex2v(8,gridCENTAUR%nHexs),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%hex2v',errorFlag)
  END IF ! errorFlag


  DO iz = 1,nz-1
    DO iy = 1,ny-1
      DO ix = 1,nx-1      
        ic = ijk2l(ix,iy,iz,nx-1,ny-1)
        
        gridCENTAUR%hex2v(1,ic) = vmap(ix+1,iy+1,iz+1)
        gridCENTAUR%hex2v(2,ic) = vmap(ix+1,iy  ,iz+1)
        gridCENTAUR%hex2v(3,ic) = vmap(ix  ,iy  ,iz+1)
        gridCENTAUR%hex2v(4,ic) = vmap(ix  ,iy+1,iz+1)
        gridCENTAUR%hex2v(5,ic) = vmap(ix+1,iy+1,iz  )
        gridCENTAUR%hex2v(6,ic) = vmap(ix+1,iy  ,iz  )
        gridCENTAUR%hex2v(7,ic) = vmap(ix  ,iy  ,iz  )
        gridCENTAUR%hex2v(8,ic) = vmap(ix  ,iy+1,iz  )        
      END DO ! ix
    END DO ! iy
  END DO ! iz

! ******************************************************************************
! Define boundaries
! ******************************************************************************

  WRITE(STDOUT,'(3X,A)') 'Determining boundary information...'

  grid3d%nBounds = 6

  ALLOCATE(grid3d%bound(grid3d%nBounds),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'grid3d%bound',errorFlag)
  END IF ! errorFlag

  grid3d%bound(1)%nVert = ny*nz
  grid3d%bound(2)%nVert = ny*nz
  grid3d%bound(3)%nVert = nx*nz
  grid3d%bound(4)%nVert = nx*nz
  grid3d%bound(5)%nVert = nx*ny
  grid3d%bound(6)%nVert = nx*ny
    
  grid3d%bound(1)%nQuads = (ny-1)*(nz-1)
  grid3d%bound(2)%nQuads = (ny-1)*(nz-1)
  grid3d%bound(3)%nQuads = (nx-1)*(nz-1)
  grid3d%bound(4)%nQuads = (nx-1)*(nz-1)
  grid3d%bound(5)%nQuads = (nx-1)*(ny-1)
  grid3d%bound(6)%nQuads = (nx-1)*(ny-1)
  
  DO ib = 1,grid3d%nBounds
    ALLOCATE(grid3d%bound(ib)%v(grid3d%bound(ib)%nVert),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN
      CALL errorHandling(ALLOCATE_ERROR,'grid3d%bound%v',errorFlag)
    END IF ! errorFlag 
    
    ALLOCATE(grid3d%bound(ib)%vMatch(grid3d%bound(ib)%nVert),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN
      CALL errorHandling(ALLOCATE_ERROR,'grid3d%bound%vMatch',errorFlag)
    END IF ! errorFlag       
  
    DO iv = 1,grid3d%bound(ib)%nVert
      grid3d%bound(ib)%vMatch(iv) = 0
    END DO ! iv
  
    ALLOCATE(grid3d%bound(ib)%quad2v(4,grid3d%bound(ib)%nQuads),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN
      CALL errorHandling(ALLOCATE_ERROR,'grid3d%bound%quad2v',errorFlag)
    END IF ! errorFlag    
  END DO ! ib

! ==============================================================================
! Boundary 1
! ==============================================================================

  WRITE(STDOUT,'(5X,A)') 'Boundary 1...'
  
  ib = 1
  ix = 1
  iq = 0
  iv = 0
  
  DO iz = 1,nz  
    DO iy = 1,ny
      iv = iv + 1
                  
      grid3d%bound(ib)%v(iv) = vmap(ix,iy,iz)  
    END DO ! iy
  END DO ! iz

  DO iz = 1,nz-1
    DO iy = 1,ny-1
      iq = iq + 1
                  
      grid3d%bound(ib)%quad2v(1,iq) = vmap(ix,iy  ,iz  )
      grid3d%bound(ib)%quad2v(2,iq) = vmap(ix,iy  ,iz+1)
      grid3d%bound(ib)%quad2v(3,iq) = vmap(ix,iy+1,iz+1) 
      grid3d%bound(ib)%quad2v(4,iq) = vmap(ix,iy+1,iz  )                        
    END DO ! iy
  END DO ! iz
  
! ==============================================================================
! Boundary 2
! ==============================================================================

  WRITE(STDOUT,'(5X,A)') 'Boundary 2...'
  
  ib = 2
  ix = nx
  iq = 0
  iv = 0
  
  DO iz = 1,nz  
    DO iy = 1,ny
      iv = iv + 1
                  
      grid3d%bound(ib)%v(iv) = vmap(ix,iy,iz)                        
    END DO ! iy
  END DO ! iz
  
  DO iz = 1,nz-1  
    DO iy = 1,ny-1
      iq = iq + 1
      
      grid3d%bound(ib)%quad2v(1,iq) = vmap(ix,iy  ,iz  )
      grid3d%bound(ib)%quad2v(2,iq) = vmap(ix,iy+1,iz  )
      grid3d%bound(ib)%quad2v(3,iq) = vmap(ix,iy+1,iz+1) 
      grid3d%bound(ib)%quad2v(4,iq) = vmap(ix,iy  ,iz+1)                        
    END DO ! iy
  END DO ! iz
  
! ==============================================================================
! Boundary 3
! ==============================================================================

  WRITE(STDOUT,'(5X,A)') 'Boundary 3...'
  
  ib = 3
  iy = 1
  iq = 0
  iv = 0
  
  DO iz = 1,nz
    DO ix = 1,nx
      iv = iv + 1
                  
      grid3d%bound(ib)%v(iv) = vmap(ix,iy,iz)                  
    END DO ! ix
  END DO ! iz

  DO iz = 1,nz-1
    DO ix = 1,nx-1
      iq = iq + 1
      
      grid3d%bound(ib)%quad2v(1,iq) = vmap(ix  ,iy,iz  )
      grid3d%bound(ib)%quad2v(2,iq) = vmap(ix+1,iy,iz  )
      grid3d%bound(ib)%quad2v(3,iq) = vmap(ix+1,iy,iz+1) 
      grid3d%bound(ib)%quad2v(4,iq) = vmap(ix  ,iy,iz+1)                        
    END DO ! ix
  END DO ! iz 
     
! ==============================================================================
! Boundary 4
! ==============================================================================

  WRITE(STDOUT,'(5X,A)') 'Boundary 4...'
  
  ib = 4
  iy = ny
  iq = 0
  iv = 0
  
  DO iz = 1,nz  
    DO ix = 1,nx
      iv = iv + 1
                  
      grid3d%bound(ib)%v(iv) = vmap(ix,iy,iz)                              
    END DO ! ix
  END DO ! iz
  
  DO iz = 1,nz-1  
    DO ix = 1,nx-1
      iq = iq + 1
      
      grid3d%bound(ib)%quad2v(1,iq) = vmap(ix  ,iy,iz  )
      grid3d%bound(ib)%quad2v(2,iq) = vmap(ix  ,iy,iz+1)
      grid3d%bound(ib)%quad2v(3,iq) = vmap(ix+1,iy,iz+1) 
      grid3d%bound(ib)%quad2v(4,iq) = vmap(ix+1,iy,iz  )                        
    END DO ! ix
  END DO ! iz
  
! ==============================================================================
! Boundary 5
! ==============================================================================

  WRITE(STDOUT,'(5X,A)') 'Boundary 5...'
  
  ib = 5
  iz = 1
  iq = 0
  iv = 0

  DO iy = 1,ny
    DO ix = 1,nx
      iv = iv + 1
                  
      grid3d%bound(ib)%v(iv) = vmap(ix,iy,iz)      
    END DO ! ix
  END DO ! iy

  DO iy = 1,ny-1
    DO ix = 1,nx-1
      iq = iq + 1
      
      grid3d%bound(ib)%quad2v(1,iq) = vmap(ix  ,iy  ,iz)
      grid3d%bound(ib)%quad2v(2,iq) = vmap(ix  ,iy+1,iz)
      grid3d%bound(ib)%quad2v(3,iq) = vmap(ix+1,iy+1,iz) 
      grid3d%bound(ib)%quad2v(4,iq) = vmap(ix+1,iy  ,iz)                        
    END DO ! ix
  END DO ! iy
  
! ==============================================================================
! Boundary 6
! ==============================================================================

  WRITE(STDOUT,'(5X,A)') 'Boundary 6...'
  
  ib = 6
  iz = nz
  iq = 0
  iv = 0
    
  DO iy = 1,ny    
    DO ix = 1,nx
      iv = iv + 1
                  
      grid3d%bound(ib)%v(iv) = vmap(ix,iy,iz)
    END DO ! ix
  END DO ! iy
  
  DO iy = 1,ny-1
    DO ix = 1,nx-1   
      iq = iq + 1
      
      grid3d%bound(ib)%quad2v(1,iq) = vmap(ix  ,iy  ,iz)
      grid3d%bound(ib)%quad2v(2,iq) = vmap(ix+1,iy  ,iz)
      grid3d%bound(ib)%quad2v(3,iq) = vmap(ix+1,iy+1,iz) 
      grid3d%bound(ib)%quad2v(4,iq) = vmap(ix  ,iy+1,iz)                        
    END DO ! ix
  END DO ! iy     

! ******************************************************************************
! Match vertices on same boundaries by check coordinate differences between 
! vertices on same patch.
! ******************************************************************************

  WRITE(STDOUT,'(3X,A)') 'Matching boundary vertices...'

  WRITE(STDOUT,'(5X,A)') 'Enter tolerance for matching:'
  READ(STDIN,*) toler

  ALLOCATE(vMatch(gridCENTAUR%nVert),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'vMatch',errorFlag)
  END IF ! errorFlag  

  nVertMatches = 0

  DO iv = 1,gridCENTAUR%nVert
    vMatch(iv) = 0
  END DO ! iv
  
! ==============================================================================
! Loop over boundaries
! ==============================================================================

  DO ib = 1,grid3d%nBounds
    nMatches = 0 
    
    WRITE(STDOUT,'(5X,A,1X,I2)') 'Boundary:',ib
    
    DO iv1 = 1,grid3d%bound(ib)%nVert
      v1 = grid3d%bound(ib)%v(iv1) 
    
      DO iv2 = iv1+1,grid3d%bound(ib)%nVert
        v2 = grid3d%bound(ib)%v(iv2) 
                
        dx = ABS(gridCENTAUR%xyz(1,v1)-gridCENTAUR%xyz(1,v2))
        dy = ABS(gridCENTAUR%xyz(2,v1)-gridCENTAUR%xyz(2,v2))
        dz = ABS(gridCENTAUR%xyz(3,v1)-gridCENTAUR%xyz(3,v2))                
                
        IF ( (dx < toler) .AND. (dy < toler) .AND. (dz < toler) ) THEN
          IF ( (grid3d%bound(ib)%vMatch(iv1) == 0) .AND. & 
               (grid3d%bound(ib)%vMatch(iv2) == 0) ) THEN                               
            nMatches = nMatches + 1
          
            grid3d%bound(ib)%vMatch(iv1) = iv2
            grid3d%bound(ib)%vMatch(iv2) = iv1
        
            IF ( vMatch(v1) == 0 .AND. vMatch(v2) == 0 ) THEN
              nVertMatches = nVertMatches + 1
                          
              vMatch(v1) = v2
              vMatch(v2) = v1
            ELSE 
              IF ( vMatch(v1) /= v2 .OR. vMatch(v2) /= v1 ) THEN 
                WRITE(*,*) 'ERROR: Inconsistent matching!'
                WRITE(*,*) v1,v2,vMatch(v1),vMatch(v2)
                STOP
              END IF ! vMatch
            END IF ! vMatch
          ELSE 
            WRITE(*,*) 'ERROR when trying to match vertices:',v1,v2
            WRITE(*,*) 'Already matched with:',grid3d%bound(ib)%vMatch(iv1), & 
                                               grid3d%bound(ib)%vMatch(iv2)
            STOP
          END IF ! grid3d%bound%vMatch
        END IF ! dx
      END DO ! iv2
    END DO ! iv 
    
    WRITE(STDOUT,'(7X,A,1X,I5)') 'Number of matches:',nMatches
  END DO ! ib

  WRITE(STDOUT,'(5X,A,1X,I5)') 'Total number of matches:',nVertMatches

  DO ib = 1,grid3d%nBounds
    DEALLOCATE(grid3d%bound(ib)%vMatch,STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN
      CALL errorHandling(DEALLOCATE_ERROR,'grid3d%bound%vMatch',errorFlag)
    END IF ! errorFlag       
  END DO ! ib
  
! ******************************************************************************
! Renumber vertices
! ******************************************************************************

  WRITE(STDOUT,'(3X,A)') 'Renumbering vertices...'
  
  ALLOCATE(vMatchInv(gridCENTAUR%nVert),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'vMatchInv',errorFlag)
  END IF ! errorFlag  
 
  DO iv = 1,gridCENTAUR%nVert ! Mark duplicated vertices
    iv2 = vMatch(iv)
    
    IF ( iv2 /= 0 ) THEN 
      IF ( iv < iv2 ) THEN 
        vMatch(iv) = 0
      ELSE 
        vMatch(iv2) = 0
      END IF ! iv      
    END IF ! vMatch
  END DO ! iv

  nVertNew = 0

  DO iv = 1,gridCENTAUR%nVert
    IF ( vMatch(iv) == 0 ) THEN 
      nVertNew = nVertNew + 1
      
      vMatch(iv) = nVertNew
      vMatchInv(nVertNew) = iv
    END IF ! vMatch  
  END DO ! iv
  
  WRITE(STDOUT,'(5X,A,I6)') 'New number of vertices:',gridCENTAUR%nVert
    
  IF ( nVertNew /= gridCENTAUR%nVert - nVertMatches ) THEN
    CALL errorHandling(NVERT_ERROR)
  END IF ! nVertNew

  gridCENTAUR%nVert = nVertNew

! ******************************************************************************
! Renumber coordinates
! ******************************************************************************

  WRITE(STDOUT,'(3X,A)') 'Renumbering coordinates...'
   
  DO iv = 1,nVertNew
    gridCENTAUR%xyz(1,iv) = gridCENTAUR%xyz(1,vMatchInv(iv))
    gridCENTAUR%xyz(2,iv) = gridCENTAUR%xyz(2,vMatchInv(iv))
    gridCENTAUR%xyz(3,iv) = gridCENTAUR%xyz(3,vMatchInv(iv))
  END DO ! iv

! ******************************************************************************
! Renumber connectivity
! ******************************************************************************

  WRITE(STDOUT,'(3X,A)') 'Renumbering connectivity...'

  DO ic = 1,gridCENTAUR%nHexs
    DO iv = 1,8
      gridCENTAUR%hex2v(iv,ic) = vMatch(gridCENTAUR%hex2v(iv,ic))
    END DO ! iv
  END DO ! iv 

! ******************************************************************************
! Renumber boundary faces
! ******************************************************************************

  WRITE(STDOUT,'(3X,A)') 'Renumbering boundary faces...'
  
  DO ib = 1,grid3d%nBounds  
    DO iq = 1,grid3d%bound(ib)%nQuads    
      DO iv = 1,4      
        v1 = grid3d%bound(ib)%quad2v(iv,iq)
                
        grid3d%bound(ib)%quad2v(iv,iq) = vMatch(v1)
      END DO ! iv
    END DO ! iq
  END DO ! ib  

  DEALLOCATE(vMatch,STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(DEALLOCATE_ERROR,'vMatch',errorFlag)
  END IF ! errorFlag 
  
  DEALLOCATE(vMatchInv,STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(DEALLOCATE_ERROR,'vMatchInv',errorFlag)
  END IF ! errorFlag   
  
! ******************************************************************************
! Remove duplicated faces by building face list for patches and checking whether
! number of original faces is equal to actual number of faces. This will differ
! if have coincident faces (obtained by renumbering faces on wake cut).  
! ******************************************************************************  

  WRITE(STDOUT,'(3X,A)') 'Removing duplicated faces from boundaries...'
  
  DO ib = 1,grid3d%nBounds  
    WRITE(STDOUT,'(5X,A,1X,I2)') 'Boundary:',ib
  
! ==============================================================================
!   Build face list. Faces which occur only once are marked with the original 
!   face index. Faces which occur twice (those along wake cut) are marked with 
!   flag -1. 
! ==============================================================================  
  
    nFaces = 0
    
    ALLOCATE(f2v(4,grid3d%bound(ib)%nQuads),STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN
      CALL errorHandling(ALLOCATE_ERROR,'f2v',errorFlag)
    END IF ! errorFlag  
      
    DO iq = 1,grid3d%bound(ib)%nQuads
      f2v(1,iq) = 0 
      f2v(2,iq) = 0    
      f2v(3,iq) = 0             
      f2v(4,iq) = 0          
    END DO ! iq  
      
    CALL CreateHashTable(grid3d%bound(ib)%nQuads)  
            
    DO iq = 1,grid3d%bound(ib)%nQuads
      fv(1) = grid3d%bound(ib)%quad2v(1,iq)
      fv(2) = grid3d%bound(ib)%quad2v(2,iq)
      fv(3) = grid3d%bound(ib)%quad2v(3,iq)
      fv(4) = grid3d%bound(ib)%quad2v(4,iq)                  
        
      CALL QuickSortInteger(fv(1:4),4)
      CALL HashBuildKey(fv(1:3),3,key)      
      CALL HashFace(iq,key,fv(1:3),nFaces,f2v,faceType)
    END DO ! iq
    
    WRITE(STDOUT,'(7X,A,I5)') 'Original number of faces:        ', &
                              grid3d%bound(ib)%nQuads
    WRITE(STDOUT,'(7X,A,I5)') 'Number of faces (w/o duplicates):',nFaces
    
    CALL DestroyHashTable

! ==============================================================================
!   If number of faces not equal to that of original grid, must generate new 
!   face list with only those faces which do not lie on wake cut.
! ==============================================================================  
    
    IF ( nFaces /= grid3d%bound(ib)%nQuads ) THEN
      cntr = 0
     
      DO iq = 1,grid3d%bound(ib)%nQuads
        iq2 = f2v(4,iq) ! Old face index
                
        IF ( iq2 > 0 ) THEN                                           
          cntr = cntr + 1      
                
          grid3d%bound(ib)%quad2v(1,cntr) = grid3d%bound(ib)%quad2v(1,iq2)
          grid3d%bound(ib)%quad2v(2,cntr) = grid3d%bound(ib)%quad2v(2,iq2)
          grid3d%bound(ib)%quad2v(3,cntr) = grid3d%bound(ib)%quad2v(3,iq2)
          grid3d%bound(ib)%quad2v(4,cntr) = grid3d%bound(ib)%quad2v(4,iq2)                        
        END IF ! iq2
      END DO ! iq
          
      grid3d%bound(ib)%nQuads = cntr
      
      WRITE(STDOUT,'(7X,A,I5)') 'Number of faces (w/o wake faces):',cntr
    END IF ! nFaces
    
    DEALLOCATE(f2v,STAT=errorFlag)
    IF ( errorFlag /= NO_ERROR ) THEN
      CALL errorHandling(DEALLOCATE_ERROR,'f2v',errorFlag)
    END IF ! errorFlag    
  END DO ! ib  

! ******************************************************************************
! Convert boundary data structure to CENTAUR format
! ******************************************************************************

  WRITE(STDOUT,'(3X,A)') 'Converting boundary data structure...'
  
  gridCENTAUR%nBounds = 6 
  
  ALLOCATE(gridCENTAUR%bInfo(3,gridCENTAUR%nBounds),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%bInfo',errorFlag)
  END IF ! errorFlag

  ALLOCATE(gridCENTAUR%bName(gridCENTAUR%nBounds),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%bName',errorFlag)
  END IF ! errorFlag  

  gridCENTAUR%nBQuads = 0 

  DO ib = 1,grid3d%nBounds
    gridCENTAUR%nBQuads = gridCENTAUR%nBQuads + grid3d%bound(ib)%nQuads
    
    gridCENTAUR%bInfo(2,ib) = 0 
    gridCENTAUR%bInfo(3,ib) = gridCENTAUR%nBQuads
  END DO ! ib  
  
  ALLOCATE(gridCENTAUR%bQuad2v(4,gridCENTAUR%nBQuads),STAT=errorFlag)
  IF ( errorFlag /= NO_ERROR ) THEN
    CALL errorHandling(ALLOCATE_ERROR,'gridCENTAUR%bQuad2v',errorFlag)
  END IF ! errorFlag  

  gridCENTAUR%nBQuads = 0
  
  DO ib = 1,grid3d%nBounds
    DO iq = 1,grid3d%bound(ib)%nQuads
      gridCENTAUR%nBQuads = gridCENTAUR%nBQuads + 1
    
      gridCENTAUR%bQuad2v(1,gridCENTAUR%nBQuads) = grid3d%bound(ib)%quad2v(1,iq)
      gridCENTAUR%bQuad2v(2,gridCENTAUR%nBQuads) = grid3d%bound(ib)%quad2v(2,iq)
      gridCENTAUR%bQuad2v(3,gridCENTAUR%nBQuads) = grid3d%bound(ib)%quad2v(3,iq)
      gridCENTAUR%bQuad2v(4,gridCENTAUR%nBQuads) = grid3d%bound(ib)%quad2v(4,iq)
    END DO ! iq
  END DO ! ib  

! ******************************************************************************
! Set patch information
! ******************************************************************************

  WRITE(STDOUT,'(3X,A)') 'Enter type of boundary for boundaries:'
  WRITE(STDOUT,'(3X,A)') '  1 - Solid wall (no-slip): Type 300'
  WRITE(STDOUT,'(3X,A)') '  2 - Solid wall (slip): Type 400'
  WRITE(STDOUT,'(3X,A)') '  3 - Symmetry: Type 500'
  WRITE(STDOUT,'(3X,A)') '  4 - Periodic: Type 600'

  DO ib = 1,grid3d%nBounds
    WRITE(STDOUT,'(5X,A,1X,I1,A)') 'Boundary ',ib,':'
    WRITE(STDOUT,'(7X,A)') 'Enter type:'
    READ(STDIN,*) gridCENTAUR%bInfo(1,ib)
    WRITE(STDOUT,'(7X,A)') 'Enter name:'
    READ(STDIN,*) gridCENTAUR%bName(ib)    
  END DO ! ib  

  WRITE(STDOUT,'(1X,A)') 'Conversion completed.'

! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE PLOT3D2CENTAUR
