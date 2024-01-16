module param_m

    implicit none

    integer, protected :: nx, nz
    integer, protected :: ntimes, ntimes_tot
    integer, protected :: t_interval ! interval between two windows
    real(8), protected :: dt
    real(8), protected :: Lx
    integer, protected :: corr_nlocs

    integer, dimension(:), allocatable, protected :: corr_loc_idx

    logical, protected :: has_input_file, uniform_step, zero_mean, corr_is_periodic, data_input_test

    character(len=128), protected :: data_filename
    character(len=128), protected :: input_filename
    character(len=128), protected :: output_spect_filename
    character(len=128), protected :: output_corr_filename

    integer, protected :: data_file
    integer, protected :: input_file
    integer, protected :: output_file

    real(8), protected :: pi

    contains

        subroutine init_param

            implicit none

            integer, parameter :: param_file = 100
            character(len=64) :: filename

            data_file = 10
            input_file = 20
            output_file = 30

            pi = acos(-1.0d0)

            write(input_filename, '(a)') 'input.dat'

            write(filename, '(a)') 'calc_spect.in'
            open(param_file, file=filename, form='formatted', action='read')

            read(param_file,*)
            read(param_file,*) nx
            read(param_file,*)
            read(param_file,*) nz
            read(param_file,*)
            read(param_file,*) ntimes
            read(param_file,*)
            read(param_file,*) ntimes_tot
            read(param_file,*)
            read(param_file,*) t_interval
            read(param_file,*)
            read(param_file,*) uniform_step
            read(param_file,*)
            read(param_file,*) dt
            read(param_file,*)
            read(param_file,*) Lx
            read(param_file,*)
            read(param_file,*) zero_mean
            read(param_file,*)
            read(param_file,*) corr_nlocs
            allocate(corr_loc_idx(corr_nlocs))
            read(param_file,*)
            read(param_file,*) corr_loc_idx(:)
            read(param_file,*)
            read(param_file,*) corr_is_periodic
            read(param_file,*)
            read(param_file,*) has_input_file
            read(param_file,*)
            read(param_file,*) data_input_test
            read(param_file,*)
            read(param_file,*) data_filename
            read(param_file,*)
            read(param_file,*) output_spect_filename
            read(param_file,*)
            read(param_file,*) output_corr_filename

            write(*,*) "nx                    = ", nx
            write(*,*) "nz                    = ", nz
            write(*,*) "ntimes                = ", ntimes
            write(*,*) "ntimes_tot            = ", ntimes_tot
            write(*,*) "t_interval            = ", t_interval
            write(*,*) "uniform_step          = ", uniform_step
            write(*,*) "dt                    = ", dt
            write(*,*) "Lx                    = ", Lx
            write(*,*) "zero_mean             = ", zero_mean
            write(*,*) "has_input_file        = ", has_input_file
            write(*,*) "data_input_test       = ", data_input_test
            write(*,*) "corr_nlocs            = ", corr_nlocs
            write(*,*) "corr_loc_idx          = ", corr_loc_idx
            write(*,*) "corr_is_periodic      = ", corr_is_periodic
            write(*,*) "data_filename         = ", trim(data_filename)
            write(*,*) "output_spect_filename = ", trim(output_spect_filename)
            write(*,*) "output_corr_filename  = ", trim(output_corr_filename)
            write(*,*) '============================='

            if (has_input_file) then
                if (.not.file_exists(input_filename)) then
                    write(*,*) "Input file does not exist!"
                end if
            end if

        end subroutine

        function file_exists(filename) result(res)

            implicit none

            character(len=*), intent(in) :: filename
            logical  :: res

            ! check if the file exists
            inquire(file=trim(filename), exist=res)

        end function

end module