module mpi_utilities

  ! Ludovico Fossa 2026
  ! performs zeroth-order interpolation of a three-dimensional cartesian field

  use config_param
  use mpi 

  implicit none

  ! MPI_COMM_WORLD variables
  integer :: ierror, rank_world, size_world
  integer :: status(MPI_STATUS_SIZE)

  ! MPI_CART variables
  integer :: rank_2d, comm2d ! comm2d is 2D communicator opaque handle in fortran
  integer :: coords(2), dims(2)
  logical :: periods(2)
  integer :: nxl_in, nyl_in, nzl_in
  integer :: nxl_out, nyl_out, nzl_out
  integer, dimension(3) :: sizes_in, sizes_out
  integer, dimension(3) :: subsizes_in, subsizes_out
  integer, dimension(3) :: starts_in, starts_out

  ! COUNTER
  integer :: i_comp

  ! I/O
  logical :: check
  integer :: field_bytes, file_bytes

  ! Subroutines
  public :: initialize_MPI, finalize_MPI, readField, writeField

contains

  subroutine initialize_MPI

    !----------------------------------------------------------------------------------------------------
    ! MPI startup (general) -----------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------
    call MPI_Init(ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank_world, ierror)
    call MPI_Comm_size(MPI_COMM_WORLD, size_world, ierror)
    ! nobody is perfect
    if((nxt_in.eq.0).or.(nyt_in.eq.0)) then 
      if(rank_world.eq.0) print*,'rank should be finite mate'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
      stop
    else if( nxt_in .ne. nyt_in ) then
      if(rank_world.eq.0) print*, 'nx and ny are different, this is above my pay grade'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
      stop
    else
      if(rank_world.eq.0) print*, &
          'reading a domain nxt_in = ', nxt_in,', nyt_in = ', nyt_in,', nzt_in = ', nzt_in
      if(rank_world.eq.0) print*, &
          'writing a domain nxt_out = ', nxt_out,', nyt_out = ', nyt_out,', nzt_out = ', nzt_out
      if(rank_world.eq.0) print*, 'running with nprocs = ', size_world 
    end if
    !----------------------------------------------------------------------------------------------------
    ! Setup 2D Cartesian topology -----------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------
    if( abs(sqrt(dble(size_world)) - int(sqrt(dble(size_world)))) .gt. 1e-10 ) then
        if(rank_world.eq.0) print*, 'the number of processors ', size_world, ' is not the square of an integer'
        if(rank_world.eq.0) print*, 'unable to generate an MPI square cartesian topology'
        call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
        stop
    else if( modulo(dble(nxt_in),sqrt(dble(size_world))).ne.0.0d0 ) then 
        if(rank_world.eq.0) print*, 'either nx = ',nxt_in,' or ',nyt_in,&
                            ' are NOT the multiple of of the number of processors ', size_world
        call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
        stop
    end if
    dims = (/ int(sqrt(dble(size_world))), int(sqrt(dble(size_world))) /)
    call MPI_Dims_create(size_world, 2, dims, ierror) ! creates 2D topology
    if (ierror /= MPI_SUCCESS) then
      print*, rank_world, "Dimension creation for Cartesian parallelization failed", ierror
      call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
    endif
    ! creates a 2d communicator --> arranges the processors on a 2d grid
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, .false., comm2d, ierror)
    if (ierror /= MPI_SUCCESS) then
      print*, rank_world, "Cartesian communicator creation failed", ierror
      call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
    endif

    ! Use the cartesian communicator from now on
    call MPI_Comm_rank(comm2d, rank_2d, ierror)
    call MPI_Cart_coords(comm2d, rank_2d, 2, coords, ierror)

    ! initialize subdomain sizes
    nxl_in = nxt_in / dims(1)
    nyl_in = nyt_in / dims(2)
    nzl_in = nzt_in
    nxl_out = nxt_out / dims(1)
    nyl_out = nyt_out / dims(2)
    nzl_out = nzt_out

    ! initialize sizes vector
    sizes_in(1) = nxt_in
    sizes_in(2) = nyt_in
    sizes_in(3) = nzt_in
    sizes_out(1) = nxt_out
    sizes_out(2) = nyt_out
    sizes_out(3) = nzt_out
    subsizes_in(1) = nxl_in
    subsizes_in(2) = nyl_in
    subsizes_in(3) = nzl_in
    subsizes_out(1) = nxl_out
    subsizes_out(2) = nyl_out
    subsizes_out(3) = nzl_out

    ! starting points in the array
    starts_in(1) = coords(1) * nxl_in
    starts_in(2) = coords(2) * nyl_in
    starts_in(3) = 0
    starts_out(1) = coords(1) * nxl_in
    starts_out(2) = coords(2) * nyl_in
    starts_out(3) = 0

    if(rank_2d.eq.0) print*,'Cartesian topology generated SUCCESSFULLY'
    if(rank_2d.eq.0) print*,'Running on MPI_CART with ',int(sqrt(dble(size_world))),' per side'

  end subroutine initialize_MPI
  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine finalize_MPI

    call MPI_Finalize(ierror)

  end subroutine finalize_MPI
  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine readField(fieldToRead,filename)

    integer :: fh ! file handle
    integer :: filetype ! file handle for the part of the file the specific processor sees
    integer(kind=MPI_OFFSET_KIND) :: displacement, total
    logical :: check
    character(len=120), intent(in) :: filename
    real(kind=8), intent(inout) :: fieldToRead(nxl_in, nyl_in, nzt_in, n_comp)
    
    if(rank_2d.eq.0) print*,'Reading from file ',trim(filename),'...'
    inquire(file=trim(filename),exist=check,size=file_bytes)
    if(check) then ! file exists
      field_bytes = int(nxt_in,MPI_OFFSET_KIND) * &
                int(nyt_in,MPI_OFFSET_KIND) * &
                int(nzt_in,MPI_OFFSET_KIND) * &
                int(n_bytes,MPI_OFFSET_KIND) ! bytes occupied by a single field
      total = field_bytes * int(n_comp,MPI_OFFSET_KIND) ! bytes occupied by a single file
      if(rank_2d.eq.0) print*, &
        'Filesize = ',file_bytes,' should not be smaller then nxt_in*nyt_in*nzt_in*n_bytes*n_comp = ', total
      if(file_bytes.ge.total) then 

        ! open file 
        call MPI_File_open( comm2d, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierror )
        call MPI_File_set_size(fh, total, ierror)   ! <-- truncates (or extends) to exactly total" bytes
        call MPI_Type_create_subarray( 3, sizes_in, subsizes_in, starts_in, &
                                      MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierror)
        call MPI_Type_commit(filetype, ierror)

        do i_comp = 1, n_comp
          ! reading component 
          displacement = int(field_bytes * (i_comp - 1), MPI_OFFSET_KIND) ! reading file from the beginning
          if(rank_2d.eq.0) print*, 'Reading field at position = ',displacement,' of ',total
          call MPI_File_set_view( fh, displacement, MPI_DOUBLE_PRECISION, filetype, "native", MPI_INFO_NULL, ierror)
          call MPI_File_read_all( fh, fieldToRead(:,:,:,i_comp), nxl_in * nyl_in * nzt_in, &
                                  MPI_DOUBLE_PRECISION, status, ierror)
        end do

        ! close file
        call MPI_File_close(fh, ierror)
        call MPI_Type_free(filetype, ierror)

      end if
    end if
    
  end subroutine readField
  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine writeField(fieldToWrite,filename)

    integer :: fh ! file handle
    integer :: filetype ! file handle for the part of the file the specific processor sees
    integer(kind=MPI_OFFSET_KIND) :: field_bytes, displacement, total
    character(len=120), intent(in) :: filename
    real(kind=8), intent(inout) :: fieldToWrite(nxl_out, nyl_out, nzt_out, n_comp)
    
    if(rank_2d.eq.0) print*,'Writing to file ',trim(filename),'...'
      
    field_bytes = int(nxt_out, MPI_OFFSET_KIND) * &
              int(nyt_out, MPI_OFFSET_KIND) * &
              int(nzt_out, MPI_OFFSET_KIND) * &
              int(n_bytes, MPI_OFFSET_KIND) ! bytes occupied by a single field
    total = field_bytes * int(n_comp, MPI_OFFSET_KIND) ! bytes occupied by a single file

    ! open file 
    call MPI_File_open( comm2d, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror )
    call MPI_File_set_size(fh, total, ierror)
    call MPI_Type_create_subarray( 3, sizes_out, subsizes_out, starts_out, &
                                  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierror)
    call MPI_Type_commit(filetype, ierror)

    do i_comp = 1, n_comp
      ! writing component 
      displacement = int(field_bytes * (i_comp - 1), MPI_OFFSET_KIND) ! reading file from the beginning
      if(rank_2d.eq.0) print*, 'Writing field at position = ', displacement,' of ', total
      call MPI_File_set_view( fh, displacement, MPI_DOUBLE_PRECISION, filetype, "native", MPI_INFO_NULL, ierror)
      call MPI_File_write_all( fh, fieldToWrite(:,:,:,i_comp), nxl_out * nyl_out * nzt_out, &
                              MPI_DOUBLE_PRECISION, status, ierror)
    end do

    ! close file
    call MPI_File_close(fh, ierror)
    call MPI_Type_free(filetype, ierror)

  end subroutine writeField

end module mpi_utilities