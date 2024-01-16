program calc_spect

  use param_m
  use data_input_m

  implicit none

  include 'fftw3.f'

  integer :: i, j, k, t, iwindow, nwindows, read_start, read_end

  real(8), allocatable, dimension(:,:) :: var, pwspt, pwspt_tot
  real(8), allocatable, dimension(:,:,:) :: data2d_window
  real(8) :: varmean
  real(8) :: varms_window_data, varms_window_spect, varms_data, varms_spect ! mean-squares

  integer(8) :: plan
  real(8), allocatable, dimension(:,:,:) :: in
  complex(8), allocatable, dimension(:,:,:) :: out

  real(8)    :: re_pw, im_pw
  character(len=64) :: output_full_filename

  integer(8) :: iret

  call init_param

  if (.not.has_input_file) then
    call open_data_files
    call user_defined_convert_input_file
    call close_data_files
  end if

  allocate(var(nx,nz))
  allocate(data2d_window(nx,ntimes,nz))
  allocate(in(nx,ntimes,nz))
  allocate(out(nx/2+1,ntimes,nz))
  allocate(pwspt(0:nx/2,(-ntimes/2+1):ntimes/2))
  allocate(pwspt_tot(0:nx/2,(-ntimes/2+1):ntimes/2))

  varmean = 0.0d0

  write(output_full_filename, '(a,a)') trim(output_spect_filename), ".txt"
  open(output_file, file=output_full_filename, form='formatted', action='write')

  ! calculate mean
  open(input_file, file=input_filename, form='unformatted', access='stream', action='read')

  write(*,*) 'reading data from', input_filename
  do j = 1, ntimes_tot
    !if (mod(j, 5000) == 0) write(*,*) 'step', j
    read(input_file) var
    varmean = varmean + sum(var)/dble(nx)/dble(nz)/dble(ntimes_tot)
  end do

  close(input_file)

  if (zero_mean) varmean = 0.0d0

  write(*,*) "mean = ", varmean

  ! read data
  open(input_file, file=input_filename, form='unformatted', access='stream', action='read')

  ! mean squares
  varms_data = 0.0d0
  varms_spect = 0.0d0

  ! calculate number of windows
  nwindows = int((ntimes_tot - ntimes)/t_interval + 1) !! one wall
  write(*,*) 'number of windows: ', nwindows
  write(*,*) '============================='

  pwspt_tot = 0.0d0

  do iwindow = 1, nwindows

    ! first window
    if (iwindow.eq.1) read_start = 1
    if (iwindow.eq.1) read_end = ntimes

    if (iwindow.gt.1) read_start = read_end + 1
    if (iwindow.gt.1) read_end = read_start + t_interval - 1

    write(*,*) 'window ', iwindow, ' new data: ', read_start, ' to ', read_end

    ! copy data of the later half of previous window to current window
    if (iwindow.gt.1) then
      data2d_window(:,1:(ntimes-t_interval),:) = data2d_window(:,(t_interval+1):ntimes,:)
    end if

    do j = 1, read_end - read_start + 1
      read(input_file) var

      if (iwindow.eq.1) data2d_window(:,j,:) = var(:,:) - varmean
      if (iwindow.gt.1) data2d_window(:,ntimes-t_interval+j,:) = var(:,:) - varmean
    end do

    varms_window_data = sum(data2d_window**2)/dble(nx)/dble(ntimes)/dble(nz)
    varms_data = varms_data + varms_window_data/dble(nwindows)


    pwspt = 0.0d0
    in = data2d_window(:,:,:)
    ! multiply by window function
    do k = 1, nz
      do t = 1, ntimes
        do i = 1, nx
          in(i,t,k) = in(i,t,k) * 0.5d0*(1.0d0-cos(2.0d0*pi*(t-1)/ntimes))
        end do
      end do
    end do

    ! thread initialization
    call dfftw_init_threads(iret)
    if (iret == 0) then
      write(*,*) 'error during thread initialization'
    end if
    call dfftw_plan_with_nthreads(16)

    call dfftw_plan_dft_r2c_3d(plan,nx,ntimes,nz,in,out,FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(plan,in,out)

    ! spect(k1, k3=0, omega)
    !do t = 1, ntimes/2+1
    !  do i = 1, nx/2+1
    !     re_pw = dreal(out(i,t,1))
    !     im_pw = dimag(out(i,t,1))
    !     pwspt(i-1,t-1) = pwspt(i-1,t-1) + (re_pw**2+im_pw**2)/dble(nx)**2/dble(ntimes)**2/dble(nz)**2
    !  end do
    !end do
    !do t = ntimes/2+2, ntimes
    !  do i = 1, nx/2+1
    !     re_pw = dreal(out(i,t,1))
    !     im_pw = dimag(out(i,t,1))
    !     pwspt(i-1,t-ntimes-1) = pwspt(i-1,t-ntimes-1) + (re_pw**2+im_pw**2)/dble(nx)**2/dble(ntimes)**2/dble(nz)**2
    !  end do
    !end do


    ! spect(k1, omega)
    do k = 1, nz
      do t = 1, ntimes/2+1
        do i = 1, nx/2+1
           re_pw = dreal(out(i,t,k))
           im_pw = dimag(out(i,t,k))
           pwspt(i-1,t-1) = pwspt(i-1,t-1) + (re_pw**2+im_pw**2)/dble(nx)**2/dble(ntimes)**2/dble(nz)**2
        end do
      end do
    end do

    do k = 1, nz
      do t = ntimes/2+2, ntimes
        do i = 1, nx/2+1
           re_pw = dreal(out(i,t,k))
           im_pw = dimag(out(i,t,k))
           pwspt(i-1,t-ntimes-1) = pwspt(i-1,t-ntimes-1) + (re_pw**2+im_pw**2)/dble(nx)**2/dble(ntimes)**2/dble(nz)**2
        end do
      end do
    end do

    varms_window_spect = 0.0d0
    do k = 1, nz
      do t = 1, ntimes
        do i = 1, nx/2+1
           re_pw = dreal(out(i,t,k))
           im_pw = dimag(out(i,t,k))
           if (i==1) then
              varms_window_spect = varms_window_spect + (re_pw**2+im_pw**2)/dble(nx)**2/dble(ntimes)**2/dble(nz)**2
           else
              varms_window_spect = varms_window_spect + 2*(re_pw**2+im_pw**2)/dble(nx)**2/dble(ntimes)**2/dble(nz)**2
           end if
        end do
      end do
    end do

    varms_spect = varms_spect + varms_window_spect/dble(nwindows)

    write(*,*) "mean-square from data (current window):", varms_window_data
    write(*,*) "mean-square from spectrum (current window):", varms_window_spect
    write(*,*) "mean-square from data:", varms_data*dble(nwindows)/dble(iwindow)
    write(*,*) "mean-square from spectrum (before renormalization):", varms_spect*dble(nwindows)/dble(iwindow)
    write(*,*) "renormalization factor:", varms_data / varms_spect
    write(*,*) '============================='

    pwspt_tot = pwspt_tot+ pwspt / dble(nwindows)
  end do


  ! renormalization
  do t = (-ntimes/2+1), ntimes/2
    do i = 0, nx/2
      !pwspt_tot(i,t) = pwspt_tot(i,t) * pw_ms_from_data / pw_ms
      pwspt_tot(i,t) = pwspt_tot(i,t) * 8.0d0/3.0d0
    end do
  end do

  write(*,*) "writing to file", output_full_filename

  do t = (-ntimes/2+1), ntimes/2
    do i = 0, nx/2
       write(output_file, '(3e24.16)') dble(i)*2.0d0*pi/Lx, dble(t)*2.0d0*pi/(ntimes*dt), pwspt_tot(i,t)*Lx*ntimes*dt/(2.0d0*pi)**2
    end do
  end do


  call dfftw_destroy_plan(plan)

  deallocate(var)
  deallocate(data2d_window)
  deallocate(in)
  deallocate(out)
  deallocate(pwspt)
  deallocate(pwspt_tot)

  close(input_file)
  close(output_file)

end program calc_spect
