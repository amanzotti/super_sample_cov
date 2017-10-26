

F90C     = gfortran -cpp
healpix = /Users/alessandromanzotti/Work/Software/Healpix_3.00
fftwdef =/Users/alessandromanzotti/Work/Software/fftw-3.3.3/api
LFLAGS  =  -llapack -lhealpix -lfftw3 -lm
FFLAGS  = -O3 -funroll-loops -ftree-vectorize -march=corei7
F90FLAGS = $(FFLAGS) -I$(healpix)/includehealpix -I$(fftwdef) -L$(cfitsio) -L$(healpix)/libhealpix $(LFLAGS) -lcfitsio 

cfitsio = /Users/alessandromanzotti/Work/Software/cfitsio




OBJFILES=  fisher_estimate.o spline.o

default: flatlens.out
all: flatlens.out

sub_flatlens.o: ./sub_flatlens.f90
	$(F90C) $(F90FLAGS) -c sub_flatlens.f90

fisher_estimate.o: ./fisher_estimate.f90 ./spline.f90
	$(F90C) $(F90FLAGS) -c flatlens.f90 sub_flatlens.f90
	
.f.o:
	f77 $(F90FLAGS) -c $<

%.o: %.f90
	$(F90C) $(F90FLAGS) -c $*.f90 

%.o: %.F90
	$(F90C) $(F90FLAGS) -c $*.F90 


flatlens.out: $(OBJFILES) 	
	$(F90C) -o flatlens.out $(OBJFILES) $(F90FLAGS)

clean:
	rm -f *.o* *.e* *.mod *.d *.pc *.obj core* *.il
