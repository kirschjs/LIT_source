F77          = gfortran

# production settings
FFLAGS  =  -O2 -mcmodel=large #-pedantic #-finteger-4-integer-8 #-funroll-loops -funroll-all-loops

# profiling setup 
#1- compile 
#2- run (without argument-file reading) 
#3- run  gprof --line ../../../src_elma/enemb.exe > tmp

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

juelma: juelma.o par/juelma
	$(F77) $(FFLAGS) -o juelma.exe juelma.o $(LIBS)

jobelma: jobelma.o par/jobelma
	$(F77) $(FFLAGS) -o jobelma.exe jobelma.o $(LIBS)

jquelma: jquelma.o par/jquelma
	$(F77) $(FFLAGS) -o jquelma.exe jquelma.o

jenelmas: jenelmas.o par/jenelmas
	$(F77) $(FFLAGS) -o jenelmas.exe jenelmas.o  $(LIBS)