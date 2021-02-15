# F77     =       g77
# FFLAGS  =       -O2 -Wall
# LIBS    =       -L/home/hmh/lib -leispack -lsu2

F77     =      gfortran

# production settings
FFLAGS  =  -Ofast -finteger-4-integer-8 -funroll-loops -funroll-all-loops -mcmodel=large
#FFLAGS  =  -O2 -finteger-4-integer-8 -mcmodel=medium
# profiling setup 
#1- compile 
#2- run (without argument-file reading) 
#3- run  gprof --line ../../../src_elma/enemb.exe > tmp
#FFLAGS  =  -Ofast -finteger-4-integer-8 -g -pg
#FFLAGS  =  -g -pg  -ffpe-summary='none' -fcray-pointer  -mcmodel=medium

LIBS    = -L. -llapack -lblas -lm -ldl

.f.o:
	$(F77) $(FFLAGS) -c -o $*.o $<


clean:
	rm -f *.o core

luelma: luelma.o par/luelma
	$(F77) $(FFLAGS) -o luelma.exe luelma.o $(LIBS)

obelma: obelma.o par/obelma
	$(F77) $(FFLAGS) -o obelma.exe obelma.o $(LIBS)

quelma: quelma.o par/quelma
	$(F77) $(FFLAGS) -o quelma.exe quelma.o 

enelmas: enelmas.o par/enelmas
	$(F77) $(FFLAGS) -o enelmas.exe enelmas.o  $(LIBS)

luise: luise.o par/LUISE
	$(F77) $(FFLAGS) -o luise.exe luise.o $(LIBS)

obem: obem.o par/OBEM
	$(F77) $(FFLAGS) -o obem.exe obem.o $(LIBS)

qual: qual.o par/QUAL
	$(F77) $(FFLAGS) -o qual.exe qual.o 

qual_M: qual_M.o par/QUAL
	$(F77) $(FFLAGS) -o qual_M.exe qual_M.o 

enemb: enemb.o par/ENEMB
	$(F77) $(FFLAGS) -o enemb.exe enemb.o $(LIBS)
