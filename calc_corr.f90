program calc_corr

  use param_m
  use data_input_m

  implicit none

  include 'mpif.h'

  integer :: ntimes_band
  integer :: i, t, nwindows, iwindow

  real(8), allocatable, dimension(:,:) :: var, p0p0, p1p1, p0p1
  real(8), allocatable, dimension(:,:,:) :: var_window
  real(8), allocatable, dimension(:) :: varmean, x_grid, t_grid, dx_grid, dt_grid

  character(len=64) :: output_full_filename

  integer :: ierr, myrank, nprocs

  call init_param

  call MPI_init(ierr)
  call MPI_comm_rank(mpi_comm_world, myrank, ierr)
  call MPI_comm_size(mpi_comm_world, nprocs, ierr)

  if (nprocs /= size(corr_loc_idx)) then
    if (myrank == 0) write(*,*) 'incorrect number of processors!'
    stop
  else
    if (myrank == 0) write(*,*) 'number of processors: ', nprocs
  end if
  write(*,*) 'rank ', myrank, ', corr_loc_idx: ', corr_loc_idx(myrank+1)
  write(*,*) '============================='

  allocate(var(nx,nz))
  allocate(var_window(nx,nz,ntimes))
  allocate(p0p0(nx,ntimes))
  allocate(p1p1(nx,ntimes))
  allocate(p0p1(nx,ntimes))
  allocate(varmean(nx))
  allocate(x_grid(nx))
  allocate(t_grid(ntimes))
  allocate(dx_grid(nx))
  allocate(dt_grid(ntimes))

  if ((.not.has_input_file) .and. (myrank == 0)) then
    call open_data_files
    call user_defined_convert_input_file
    call close_data_files
  end if

  call MPI_barrier(mpi_comm_world, ierr)

  if (myrank == 0) then
    call user_defined_grid_file(x_grid, t_grid)
  end if

  call MPI_bcast(x_grid, nx, MPI_REAL8, 0, mpi_comm_world, ierr)
  call MPI_bcast(t_grid, ntimes, MPI_REAL8, 0, mpi_comm_world, ierr)

  ntimes_band = ntimes / 2 ! ntimes should be an even number

  varmean = 0.0d0

  ! calculate mean
  if ((myrank == 0).and.(.not.zero_mean)) then
    open(input_file, file=input_filename, form='unformatted', access='stream', action='read')
    write(*,*) 'reading data from', input_filename

    do t = 1, ntimes_tot
      if (mod(t,5000) == 0) write(*,*) 'step', t
      read(input_file) var
      varmean = varmean + sum(var,dim=2)/dble(nz)/dble(ntimes_tot)
    end do

    close(input_file)
  end if

  if (myrank == 0) then
    write(*,*) "mean = ", sum(varmean)/dble(nx)
    write(*,*) '============================='
  end if

  call MPI_bcast(varmean, nx, MPI_REAL8, 0, mpi_comm_world, ierr)

  ! calculate number of windows
  if (t_interval > 0) then
    nwindows = int((ntimes_tot - ntimes)/t_interval+1)
  else
    nwindows = 1
  end if

  if (myrank == 0) then
    write(*,*) 'number of windows: ', nwindows
    write(*,*) '============================='
  end if

  open(input_file, file=input_filename, form='unformatted', access='stream', action='read')

  p0p0 = 0.0d0
  p1p1 = 0.0d0
  p0p1 = 0.0d0

  do iwindow = 1, nwindows

    if (myrank == 0) write(*,*) 'window', iwindow, ':', 1+t_interval*(iwindow-1), 'to', ntimes+t_interval*(iwindow-1)

    if (iwindow == 1) then
      do t = 1, ntimes
        read(input_file) var_window(:,:,t)
      end do
    else
      ! copy data from the previous window
      do t = 1, ntimes-t_interval
        var_window(:,:,t) = var_window(:,:,t+t_interval)
      end do
      ! read new data
      do t = 1, t_interval
        read(input_file) var_window(:,:,ntimes-t_interval+t)
      end do
    end if

    do t = 1, ntimes
      do i = 1, nx
        p0p0(i,t) = p0p0(i,t) + sum((var_window(corr_loc_idx(myrank+1),:,ntimes_band) - varmean(corr_loc_idx(myrank+1)))**2) / dble(nz) /dble(nwindows)
        p1p1(i,t) = p1p1(i,t) + sum((var_window(i,:,t) - varmean(i))**2) / dble(nz) / dble(nwindows)
        p0p1(i,t) = p0p1(i,t) + sum((var_window(corr_loc_idx(myrank+1),:,ntimes_band) - varmean(corr_loc_idx(myrank+1)))*(var_window(i,:,t) - varmean(i))) / dble(nz) /dble(nwindows)
      end do
    end do
  end do

  close(input_file)

  do i = 1, nx
    dx_grid(i) = x_grid(i) - x_grid(corr_loc_idx(myrank+1))

    !for periodic domain
    if (corr_is_periodic) then
      if (i-corr_loc_idx(myrank+1) > nx/2) then
        dx_grid(i) = dx_grid(i) - Lx
      else if (i-corr_loc_idx(myrank+1) < -nx/2) then
        dx_grid(i) = dx_grid(i) + Lx
      end if
    end if
  end do

  do t = 1, ntimes
    if (uniform_step) then
      dt_grid(t) = (t - ntimes_band) * dt
    else
      dt_grid(t) = t_grid(t) - t_grid(ntimes_band)
    end if
  end do

  if (myrank == 0) write(*,*) '============================='
  write(output_full_filename, '(a,a,i4.4,a)') trim(output_corr_filename), "_", corr_loc_idx(myrank+1), ".txt"
  write(*,*) "writing to file ", output_full_filename
  open(output_file, file=output_full_filename, form='formatted', action='write')
  do t = 1, ntimes
    do i = 1, nx
      write(output_file, '(3e24.16)') dx_grid(i), dt_grid(t), p0p1(i,t)/(sqrt(p0p0(i,t))*sqrt(p1p1(i,t)))
      !write(output_file, '(3e24.16)') dx_grid(i), dt_grid(t), p0p1(i,t)
    end do
  end do
  close(output_file)

  call MPI_finalize(ierr)

  deallocate(var)
  deallocate(var_window)
  deallocate(p0p0)
  deallocate(p1p1)
  deallocate(p0p1)
  deallocate(varmean)
  deallocate(x_grid)
  deallocate(t_grid)
  deallocate(dx_grid)
  deallocate(dt_grid)

end program calc_corr
