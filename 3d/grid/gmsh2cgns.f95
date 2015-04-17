!----------------------------------------------------------------------
!! Converts gmsh file into cgns file for use with freecfd
!! works only for tetrahedra with triangular boundary faces
!! Compilation:
!!    gfortran gmsh2cgns.f95 -I/usr/local/include -L/usr/local/lib -lcgns
!! Written by Praveen. C, <http://math.tifrbng.res.in/~praveen>
!----------------------------------------------------------------------
PROGRAM GMSH_TO_CGNS

  !--------------------------------------------------
  IMPLICIT NONE

  include 'cgnslib_f.h'  
  !--------------------------------------------------
  !! Number of vertices and associated index
  INTEGER :: ns, is

  !! Number of elements (faces + tetrahedra) and associated idex
  INTEGER :: nelem, iel

  !! Number of tetrahedra and associated index
  INTEGER :: nt, jt

  !! Number of faces and associated index
  INTEGER :: nfac, ifac

  !! Type of gmsh element and array containing
  !! the number of nodes belonging to this 
  !! element according its type.
  INTEGER :: type_elem
  INTEGER, DIMENSION(4) :: node_per_element
  
  !! Type of physical group which are used to
  !! determine the type of boundary conditions
  !! and associated index
  INTEGER :: type_phys, itype

  !! Number of boundary logicals and associated index
  INTEGER :: nitab, itab

  !! Name of the input and output files
  CHARACTER(LEN=256) :: fname_mesh
  CHARACTER(LEN=256) :: fname_out

  !! String used to skip some parts of the input file
  CHARACTER(LEN=256) :: string

  !! Array of node coordinates
  REAL(KIND=KIND(1.D0)), DIMENSION(:), POINTER :: x, y, z

  !! Array of connectivity tetrahedra -> nodes
  !! The first index denotes the id of the nodes 
  !! The second one denotes the id of the tetrahedron
  INTEGER, DIMENSION(:,:), POINTER :: nu

  !! Array of connectivity faces -> nodes
  !! The first index denotes the id of the nodes
  !! The second one denotes the id of the face
  INTEGER, DIMENSION(:,:), POINTER :: nsfac

  !! Array containing the boundary logicals for each face
  INTEGER, DIMENSION(:), POINTER :: logfac

  !! Array of correspondance between
  !! the gmsh physical type and the boundary type
  !! First index denotes the current boundary
  !! Second index denotes:
  !!   -the gmsh type if it is equal to 1
  !!   -the numesis type if it is equal to 2
  INTEGER, DIMENSION(:,:), POINTER :: Tab
  
  !! Provisional array storing the id of the nodes
  INTEGER :: list(4)

  !! loop index
  INTEGER :: i

  !! Integer which denotes a useless data read in the input file
  INTEGER :: inull
  INTEGER :: ntags
  
  integer :: isize(3,3), ipnts(2)
  character(len=32) :: basename,zonename
  integer :: icelldim, iphysdim, index_file, index_base, index_coord, ier
  integer :: index_zone, index_section, index_bc
  integer :: ielem_no, nelem_start, nelem_end, nbdyelem, icount
  integer, allocatable :: felem(:,:)
  character(len=16) :: bcname

  !-------------------------------------------------
  ! Get arguments of the command line
  ! 1- input file
  ! 2- output file
  CALL GETARG(1,fname_mesh)

  ! Verification of the command line
  IF(TRIM(fname_mesh) == "")THEN     
     PRINT '(" ")'
     PRINT '("This convertor requires one argument")'
     PRINT '(" >>> ./gmsh2cgns  filename")'
     PRINT '(" ")'
     STOP
  END IF


  !-----------------
  ! read the mesh 
  !-----------------
  PRINT '(" >>> GMSH file reading ... ")'

  OPEN(UNIT=10, FILE=TRIM(fname_mesh)//'.msh', STATUS='OLD')

  ! Skip some lines
  string = " "
  DO WHILE (TRIM(string) /= "$Nodes")
     READ(UNIT = 10,  FMT = *) string
  END DO

  ! number of nodes and allocation of the array
  READ(UNIT = 10, FMT = *) ns
  
  ALLOCATE(x(ns))
  ALLOCATE(y(ns))
  ALLOCATE(z(ns))

  ! read the node coordinates
  DO is = 1, ns
     READ(UNIT = 10, FMT = *) inull, x(is), y(is), z(is)
  END DO

  ! Skip some lines
  string = " "
  DO WHILE (TRIM(string) /= "$Elements")
     READ(UNIT = 10, FMT = *) string
  END DO

  ! number of faces and elements and allocation of the arrays
  READ(UNIT = 10, FMT = *) nelem

  ALLOCATE(nu(4,nelem))
  ALLOCATE(nsfac(3,nelem))
  ALLOCATE(logfac(nelem))

  ! Correspondance between type of element and number of nodes
  ! 1 -> unknown
  ! 2 -> triangle
  ! 3 -> unknown
  ! 4 -> tetrahedron
  node_per_element(1)=0
  node_per_element(2)=3
  node_per_element(3)=0
  node_per_element(4)=4

  ! read elements and storage according their types
  nt=0
  nfac =0

  DO iel=1,nelem

     READ(UNIT   = 10, &
          FMT    = *) inull, type_elem, ntags, type_phys, inull, inull, &
          (list(i), i=1,node_per_element(type_elem))

     IF(ntags /= 3) THEN
        PRINT*,'!!! ntags is not equal to 3'
        PRINT*,'Element =',iel
        PRINT*,'ntags   =',ntags
        STOP
     ENDIF

     ! case : tetrahedron
     IF(type_elem == 4) THEN

        nt = nt + 1

        nu(1,nt) = list(1)
        nu(2,nt) = list(2)
        nu(3,nt) = list(3)
        nu(4,nt) = list(4)

        ! case : triangle
     ELSE IF(type_elem == 2) THEN

        nfac = nfac + 1

        logfac(nfac) = type_phys
        nsfac(1,nfac) = list(1)
        nsfac(2,nfac) = list(2)
        nsfac(3,nfac) = list(3)

        ! else
     ELSE
        PRINT '("bad type of element")'
        STOP
     END IF

  END DO

  CLOSE(UNIT = 10)

  PRINT '(" ")'
  PRINT '(" >>>  End of reading ")'
  PRINT '(" ")'
  PRINT '(" >>>  Nb of elements: ", i8)', nt
  PRINT '(" >>>  Nb of faces: ", i8)', nfac

  !--------------------------
  ! treatment of boundaries
  !--------------------------
  ALLOCATE(Tab(1000,4))
  
  ! modify the logic ID
  nitab=0
  
  DO ifac=1,nfac
  
     type_phys = logfac(ifac)
  
     ! search in the table of boundary types
     itype=0
     DO itab=1,nitab
        IF(Tab(itab,1) == type_phys)THEN
          itype=1
          Tab(itab,2) = Tab(itab,2) + 1
        ENDIF
     END DO
  
     ! if not found : creation of a new type
     IF(itype == 0) THEN
  
        nitab = nitab+1
        Tab(nitab,1) = type_phys
        Tab(nitab,2) = 1
  
     END IF
  
  END DO
  
  PRINT '(" >>>  Nb of boundary face types: ", i8)', nitab
  DO i=1,nitab
    PRINT*, Tab(i,1), Tab(i,2)
  ENDDO

  !---------------------
  ! write CGNS file 
  !---------------------
  fname_out = TRIM(fname_mesh)//'.cgns'
  PRINT*,'     CGNS file written into ', TRIM(fname_out)
  ! open CGNS file for write
      call cg_open_f(fname_out,MODE_WRITE,index_file,ier)
  ! create base (user can give any name)
      basename='Base'
      icelldim=3
      iphysdim=3
      call cg_base_write_f(index_file,basename,icelldim,iphysdim, &
                           index_base,ier)
  ! define zone name (user can give any name)
      zonename = 'Zone1'
  ! no of vertex
      isize(1,1)=ns
  ! no of cell
      isize(2,1)=nt
  ! boundary vertex size (zero if elements not sorted)
      isize(3,1)=0
  ! create zone
      call cg_zone_write_f(index_file,index_base,zonename,isize, &
                           Unstructured,index_zone,ier)
  ! write grid coordinates (user must use SIDS-standard names here)
      call cg_coord_write_f(index_file,index_base,index_zone,RealDouble, &
                            'CoordinateX',x,index_coord,ier)
      call cg_coord_write_f(index_file,index_base,index_zone,RealDouble, &
                            'CoordinateY',y,index_coord,ier)
      call cg_coord_write_f(index_file,index_base,index_zone,RealDouble, &
                            'CoordinateZ',z,index_coord,ier)
  ! set element connectivity:
  ! do all the HEX_8 elements (this part is mandatory):
  ! maintain SIDS-standard ordering
      ielem_no=0
  ! index no of first element
      nelem_start=1
  ! index no of last element
      nelem_end=nt
  ! unsorted boundary elements
      nbdyelem=0
  ! write tetrahedra element connectivity (user can give any name)
      call cg_section_write_f(index_file,index_base,index_zone, &
                              'TET',TETRA_4,nelem_start,nelem_end,nbdyelem,nu, &
                              index_section,ier)
  ! add boundary elements
  DO i=1,nitab
    icount = 0
    allocate(felem(3,Tab(i,2)))
    nelem_start = nelem_end + 1
    nelem_end   = nelem_end + Tab(i,2)
    Tab(i,3) = nelem_start
    Tab(i,4) = nelem_end
    DO ifac=1,nfac
      IF(logfac(ifac)==Tab(i,1))THEN
        icount = icount + 1
        felem(1,icount) = nsfac(1,ifac)
        felem(2,icount) = nsfac(2,ifac)
        felem(3,icount) = nsfac(3,ifac)
      ENDIF
    ENDDO
    call setbcname(Tab(i,1), bcname)
    call cg_section_write_f(index_file,index_base,index_zone, &
                            bcname,TRI_3,nelem_start,nelem_end,nbdyelem,felem, &
                            index_section,ier)
    deallocate(felem)
  ENDDO
  ! add boundary conditions
  DO i=1,nitab
    icount = 2
    ipnts(1) = Tab(i,3)
    ipnts(2) = Tab(i,4)
    call setbcname(Tab(i,1), bcname)
    print '(" >>>  BC name = ",a8)', bcname
    call cg_boco_write_f(index_file,index_base,index_zone,bcname, &
                         Null,ElementRange,icount,ipnts,index_bc,ier)
  ENDDO
  ! close CGNS file
      call cg_close_f(index_file,ier)

  ! Deallocation
  DEALLOCATE(x)
  DEALLOCATE(y)
  DEALLOCATE(z)

  DEALLOCATE(nu)
  DEALLOCATE(logfac)
  DEALLOCATE(nsfac)
  DEALLOCATE(Tab)

  PRINT '(" ")'
  PRINT '(" >>>  Finished ! ")'

END PROGRAM GMSH_TO_CGNS

!------------------------------------------------------------------------------
subroutine setbcname(phytag, bcname)
  implicit none
  integer :: phytag
  character(len=*) :: bcname
  
  if(phytag < 10)then
    write(unit=bcname,fmt='(i1)') phytag
  elseif(phytag < 100)then
    write(unit=bcname,fmt='(i2)') phytag
  elseif(phytag < 1000)then
    write(unit=bcname,fmt='(i3)') phytag
  elseif(phytag < 10000)then
    write(unit=bcname,fmt='(i4)') phytag
  elseif(phytag < 100000)then
    write(unit=bcname,fmt='(i5)') phytag
  elseif(phytag < 1000000)then
    write(unit=bcname,fmt='(i6)') phytag
  else
    print*,'Error in setbcname'
    stop
  endif
  
end subroutine setbcname
