program interpolate_MPI

  ! Ludovico Fossa 2026
  ! performs zeroth-order interpolation of a three-dimensional cartesian field
  use config_param
  use mpi_utilities

  implicit none

  ! VARIABLES 
  integer :: aspect_ratio
  character(len=20) :: iStepChar
  character(len=120) :: filename, vtkname
  real(kind=8), allocatable, dimension(:,:,:,:) :: field_in, field_out
  integer :: i_in, j_in, k_in, i_out, j_out, k_out, i_stride, j_stride, k_stride

  ! START PROGRAM
  call initialize_MPI
  if (rank_2d.eq.0) write(6,*) 'interpolate_MPI.f90 start'

  ! CHECK WHETHER ASPECT RATIO IS INTEGER AND CONSISTENT
  if((mod(nxt_out,nxt_in).gt.0).or.(mod(nyt_out,nyt_in).gt.0).or.(mod(nzt_out,nzt_in).gt.0)) then 
    if(rank_2d.eq.0) print*, 'aspect ratio is not an integer'
    call mpiEnd(.true.)
    stop 
  else if((nxt_in * nyt_out).eq.(nxt_out * nyt_in) .and. &
          (nxt_in * nzt_out).eq.(nxt_out * nzt_in) .and. &
          (nzt_in * nyt_out).eq.(nzt_out * nyt_in)) then
    aspect_ratio = nxt_out / nxt_in
    if(rank_2d.eq.0) print*, 'aspect ratio = ',aspect_ratio
  endif

  !-------------------------------------------
  ! READ AND INTERPOLATE FIELDS
  !-------------------------------------------
  ! ALLOCATE
  allocate(field_in(nxl_in,nyl_in,nzt_in,n_comp))
  write(iStepChar,'(i12.12)') iStep
  filename = trim(datadir)//trim(prefix)//trim(iStepChar)//trim(suffix)
  if(rank_2d.eq.0) write(6,*) 'reading ',trim(filename)
  inquire(file=trim(filename),exist=check,size=file_bytes)
  if(rank_2d.eq.0) write(6,*) 'exist ',check
  if(check) then ! file exists
    call readField(field_in,filename)
  else
    if(rank_2d.eq.0) write(*,*) 'Missing/wrong input file ',trim(filename)
    !call mpiEnd(.FALSE.)
    if(rank_2d.eq.0) write(*,*) 'stop'
    deallocate(field_in)
    stop
  endif

  vtkname = 'snapshot_in'
  call printVTK_slice(field_in(:,:,:,1),vtkname,nxl_in,nyl_in,nzl_in)

  nxl_out = nxl_in * aspect_ratio
  nyl_out = nyl_in * aspect_ratio
  nzl_out = nzt_in * aspect_ratio
  allocate(field_out(nxl_out,nyl_out,nzt_out,n_comp))
  
  do i_comp = 1, n_comp
    do k_in = 1, nzt_in
      k_out = (k_in - 1) * 2 + 1
      do j_in = 1, nyl_in
        j_out = (j_in - 1) * 2 + 1
        do i_in = 1, nxl_in
          i_out = (i_in - 1) * 2 + 1
          if(aspect_ratio.gt.1) then ! interpolate on a larger field
            do i_stride = 0, aspect_ratio - 1
              do j_stride = 0, aspect_ratio - 1
                do k_stride = 0, aspect_ratio - 1
                  field_out(i_out+i_stride,j_out+j_stride,k_out+k_stride,i_comp) = & 
                            field_in(i_in,j_in,k_in,i_comp)
                end do
              end do
            end do
          end if
        end do
      end do
    end do 
  end do
  deallocate(field_in)

  !-------------------------------------------
  ! WRITE
  !-------------------------------------------
  filename = trim(datadir)//trim(prefix)//trim(iStepChar)
  if(rank_2d.eq.0) print*, 'Writing on ', trim(filename),'...'
  call writeField(field_out,filename)

  vtkname = 'snapshot_out'
  call printVTK_slice(field_out(:,:,:,1),vtkname,nxl_out,nyl_out,nzl_out)

  deallocate(field_out)
  
  call finalize_MPI

end program interpolate_MPI
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine printVTK_slice(domain_local,outputname,size_x,size_y,size_z)

  ! print 2D slice in vtk (verify that the output is consistent and correct)
  use config_param
  use mpi_utilities

  implicit none

  integer :: i, j, slice, istart, jstart, iend, jend
  character(len=120) :: filename, outputname
  character(len=12) :: istepChar
  integer, intent(in) :: size_x, size_y, size_z
  real(kind=8), intent(in) :: domain_local(size_x,size_y,size_z)
  real(kind=8), allocatable, dimension(:,:) :: domain_global
  real(kind=8), allocatable, dimension(:,:) :: buf
  integer :: src, src_coords(2)
  integer :: size_x_tot, size_y_tot

  slice = int(size_z/2)
  if(rank_2d.eq.0) print*,'printing slice = ',slice,'/',size_z
  if(rank_2d.eq.0) print*,'checking size world = ',size_world
  if(rank_2d.eq.0) print*,'outputname = ',outputname

  size_x_tot = size_x * dims(1)
  size_y_tot = size_y * dims(2)

  ! allocate receive buffer for root
  allocate(domain_global(size_x_tot,size_y_tot))
  domain_global = 0.0
  
  ! allocate 2d buffer (slice) for sending 2D slice
  allocate(buf(size_x,size_y))
  buf = domain_local(:,:,slice)

  ! send/recv to root
  if (rank_2d.eq.0) then

    ! roots copies own block
    istart = 1 + coords(1) * size_x
    iend   = istart + size_x - 1
    jstart = 1 + coords(2) * size_y
    jend   = jstart + size_y - 1
    domain_global(istart:iend,jstart:jend) = buf

    ! receive from all other ranks
    do src = 0, size_world-1
      
      if (src == 0) cycle

      call MPI_Recv(src_coords, 2, MPI_INTEGER, src, 100, comm2d, status, ierror)

      istart = 1 + src_coords(1) * size_x
      iend   = istart + size_x - 1
      jstart = 1 + src_coords(2) * size_y
      jend   = jstart + size_y - 1

      call MPI_Recv(domain_global(istart:iend, jstart:jend), size_x * size_y, MPI_DOUBLE_PRECISION, &
                        src, 101, comm2d, status, ierror)

    end do
  
  else

    ! send coords first, data later
    call MPI_Send(coords, 2, MPI_INTEGER, 0, 100, comm2d, ierror)
    call MPI_Send(buf, size_x * size_y, MPI_DOUBLE_PRECISION, 0, 101, comm2d, ierror)
  
  end if

  deallocate(buf)

  ! writing the final vtk file here and deallocating domain_global
  if(rank_2d.eq.0) then
    write(iStepChar,'(i12.12)') iStep
    filename=trim(datadir)//trim(prefix)//'_'//trim(outputname)//'_'//iStepChar//'.vtk'
    open(unit=10, file=filename, status='replace')

    ! --- VTK header ---
    write(10,'(A)') "# vtk DataFile Version 3.0"
    write(10,'(A)') "Enstrophy test grid"
    write(10,'(A)') "ASCII"
    write(10,'(A)') "DATASET STRUCTURED_POINTS"
    write(10,'(A,3I6)') "DIMENSIONS", size_x_tot, size_y_tot, 1
    write(10,'(A,3F12.5)') "ORIGIN", 0.0, 0.0, 0.0
    write(10,'(A,3F12.5)') "SPACING", 1.0, 1.0, 1.0
    write(10,'(A,I12)') "POINT_DATA", size_x_tot * size_y_tot * 1
    write(10,'(A)') "SCALARS scalars double"
    write(10,'(A)') "LOOKUP_TABLE default"

    ! --- Write values ---
    do j=1,size_y_tot
      do i=1,size_x_tot
        write(10,'(E16.8)') domain_global(i,j)
      end do
    end do

    close(10)
    print *, "Wrote ", trim(filename)

  end if

  deallocate(domain_global)

end subroutine printVTK_slice