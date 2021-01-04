F77     =      mpif77
CC      =      mpicc

FFLAGS  =       -O2 -cpp -mcmodel=medium
#FFLAGS  =       -Ofast -cpp -funroll-loops -funroll-all-loops -mcmodel=medium  

# -g -pg  -ffpe-summary='none' -fcray-pointer  -mcmodel=medium
#-fno-align-commons#-fdefault-integer-8 -ffpe-trap=denormal  -fno-align-commons 

LIBS    = -L. -llapack -lblas

#FFLAGS  =       -Ofast -g -traceback -i8 -mcmodel=medium
#FFLAGS =  -ipo -O3 -no-prec-div -fp-model fast=2 -xHost # = fast w/o -static
#FFLAGS  =      -O2 -g -traceback -i8 
#LIBS    = -L${MKLROOT} -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lm -ldl

.f.o:
	$(F77) $(FFLAGS) -c -o $*.o $<

clean:
	rm -f *.o core

qual: qual.o
	$(F77) $(FFLAGS) -o qual.exe qual.o $(LIBS)

luise: luise.o
	$(F77) $(FFLAGS) -o luise.exe luise.o $(LIBS)

obem: obem.o
	$(F77) $(FFLAGS) -o obem.exe obem.o $(LIBS)

qual_joki: qual_joki.o
	$(F77) $(FFLAGS) -o qual_joki.exe qual_joki.o $(LIBS)

qual_joki_2: qual_joki_2.o
	$(F77) $(FFLAGS) -o qual_joki_2.exe qual_joki_2.o $(LIBS)

luise_joki: luise_joki.o
	$(F77) $(FFLAGS) -o luise_joki.exe luise_joki.o $(LIBS)

obem_joki: obem_joki.o
	$(F77) $(FFLAGS) -o obem_joki.exe obem_joki.o $(LIBS)

enemb: enemb.o
	$(F77) $(FFLAGS) -o enemb.exe enemb.o $(LIBS)

tdenemb: tdenemb.o
	$(F77) $(FFLAGS) -o tdenemb.exe tdenemb.o $(LIBS)

sammelmag: sammelmag.o
	$(F77) $(FFLAGS) -o sammelmag.exe sammelmag.o $(LIBS)	
