all: tn.a

GOBJS =	blas.o tn.o
FFLAGS = -g -c -FI -zero -fno-automatic -std=legacy -fPIC
FC = gfortran

.SUFFIXES: .f .o

tn.a: $(GOBJS)
	ar rs $@ $+

clean:
	rm -f -- *.o *.a
