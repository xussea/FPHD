FC=gfortran 
FFLAGS=-O2 -fopenmp -fdefault-real-8
LAPACK=-llapack -lblas
OUT=FPHD

all:principal clean

principal: subroutines
	${FC} ${FFLAGS} FPHD_main.f90 FPHD_readsub.o FPHD_sttanalysis.o FPHD_exactcharges.o FPHD_tdmjoinjacgetin.o ${LAPACK} -o ${OUT}
	mkdir -p ../bin
	mv ${OUT} ../bin/
subroutines:
	${FC} ${FFLAGS} -c FPHD_readsub.f90 FPHD_sttanalysis.f90 FPHD_exactcharges.f90 FPHD_tdmjoinjacgetin.f90
clean:
	-rm -f *.o *.mod

