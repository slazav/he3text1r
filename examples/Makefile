all: test_c test_f77 test_f90

CFLAGS=-Wall
LIBS= ../libtext1r.a ../../he3lib/libhe3.a ../tn/tn.a
FC=gfortran

test_c:  test_c.o $(LIBS)
	$(FC) -g $+ -lgfortran -lm -o $@

test_f77: test_f77.o $(LIBS)
	$(FC) -g $+ -o $@

test_f90: test_f90.o $(LIBS)
	$(FC) -g $+ -o $@


%.o: %.f90
	$(FC) -g $(FFLAGS) -c $< -o $@

clean:
	rm -f test_c test_f77 test_f90 *.o
