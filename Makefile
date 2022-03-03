include ./compiler.bridges2.inc

#all: dgemv

dgemv: dgemv.o
	$(FC) -o $@ $^ ${LDFLAGS}

%.o : %.f
	$(FF) $(FFLAGS) -I${MKLINCLUDE} -I${MKLINCLUDE}/intel64/lp64 -c $<

clean:
	rm dgemv.o dgemv
	

