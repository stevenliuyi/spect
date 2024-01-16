module data_input_m

    use param_m

    implicit none

    contains

        subroutine open_data_files

            implicit none

            open(data_file, file=data_filename, form='unformatted', access='stream', action='read')

            open(input_file, file=input_filename, form='unformatted', access='stream', action='write')

        end subroutine


        subroutine close_data_files

            implicit none

            close(data_file)
            close(input_file)

            write(*,*) 'input data saved at', input_filename
            write(*,*) '============================='

        end subroutine

        subroutine user_defined_convert_input_file

            implicit none

            integer, parameter :: ny = 192
            real(8), allocatable, dimension(:,:) :: var, var_uniform
            integer :: i, j, k, idx

            character(len=64) :: grid_filename
            integer :: grid_file, nblocks
            integer, dimension(3) :: size_tot
            real(8), allocatable, dimension(:,:,:,:) :: xyz
            real(8), allocatable, dimension(:) :: x_grid, x_uniform
            integer, allocatable, dimension(:) :: idx_uniform

            allocate(var(nx+1,nz+1))
            !allocate(var_uniform(nx+1,nz+1))
            allocate(xyz(nx+1,ny,nz+1,3))
            allocate(x_grid(nx+1))
            !allocate(x_uniform(nx))
            !allocate(idx_uniform(nx))

            grid_file = 50
            write(grid_filename, '(a)') 'grid.3d.dat'
            open(grid_file, file=grid_filename, form='unformatted', access='stream', action='read')

            ! PLOT3D format
            read(grid_file) nblocks
            read(grid_file) size_tot
            read(grid_file) xyz(:,:,:,1)
            read(grid_file) xyz(:,:,:,2)
            read(grid_file) xyz(:,:,:,3)

            x_grid = xyz(:,ny,1,1)
            deallocate(xyz)

            close(grid_file)

            ! directly save without modification
            !do j = 1, ntimes_tot
            !    read(data_file) var
            !end do
            do j = 1, ntimes_tot
                read(data_file) var
                write(input_file) var(1:nx,1:nz)
            end do

            ! interpolation
            !do i = 1, nx+1
            !    x_uniform(i) = x_grid(1) + (x_grid(nx+1) - x_grid(1)) / (nx) * (i-1)
            !end do

            !do i = 1, nx+1
            !    idx_uniform(i) = minloc(x_uniform(i) - x_grid, dim=1, mask=(x_uniform(i) >= x_grid))
            !end do
            !idx_uniform = max(min(idx_uniform, nx), 1)

            !do j = 1, ntimes_tot
            !    write(*,*) 'step', j
            !    read(data_file) var

            !    do k = 1, nz+1
            !        do i = 1, nx+1
            !            idx = idx_uniform(i)
            !            var_uniform(i,k) = var(idx,k) + (x_uniform(i) - x_grid(idx)) * (var(idx+1,k) - var(idx,k)) / (x_grid(idx+1) - x_grid(idx))
            !        end do
            !    end do

            !    write(input_file) var_uniform(1:nx,1:nz)
            !end do

            deallocate(var)
            !deallocate(var_uniform)
            deallocate(x_grid)
            !deallocate(x_uniform)
            !deallocate(idx_uniform)

        end subroutine

        subroutine user_defined_grid_file(x_grid, t_grid)

            implicit none

            real(8), intent(out) :: x_grid(1:nx), t_grid(1:ntimes)

            integer, parameter :: ny = 192

            character(len=64) :: grid_filename

            integer :: grid_file, nblocks
            integer, dimension(3) :: size_tot
            real(8), allocatable, dimension(:,:,:) :: xyz

            allocate(xyz(nx,ny,nz))

            grid_file = 50
            write(grid_filename, '(a)') 'grid.3d.dat'
            open(grid_file, file=grid_filename, form='unformatted', action='read')

            ! PLOT3D format
            read(grid_file) nblocks
            read(grid_file) size_tot
            read(grid_file) xyz(:,:,:)

            x_grid = xyz(:,ny,1)

            close(grid_file)
            deallocate(xyz)

            t_grid = 0.0d0

        end subroutine


end module
