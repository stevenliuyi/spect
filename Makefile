F90 = mpif90

INCLUDE = -I $(CSE_FFTW_INCLUDE_DIR)
LIBS = -L $(CSE_FFTW_LIB_DIR) -lfftw3 -lfftw3_threads -lm -lpthread

OPTIONS = -ffree-line-length-none
#OPTIONS = -ffree-line-length-none -Wall
#OPTIONS = -mcmodel=medium -heap-arrays -warn all,noexternal

OPTIM = -O3

F90FLAGS = $(OPTIM)
LDFLAGS = $(OPTIM)


OBJ = param.o data_input.o
SPECT_OBJ = $(OBJ) calc_spect.o
CORR_OBJ = $(OBJ) calc_corr.o

all: calc_spect calc_corr

calc_spect: $(SPECT_OBJ)
	$(F90) -o $@ $(SPECT_OBJ) $(LIBS)

calc_corr: $(CORR_OBJ)
	$(F90) -o $@ $(CORR_OBJ) $(LIBS) -lmpi

clean:
	rm -f *.o *.mod *~ calc_spect

%.o : %.f90
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) -c $<
