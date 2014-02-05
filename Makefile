FFLAGS= -Werror -Wconversion\
  -Wline-truncation\
  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\
  -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow\
  -Wno-unused-parameter -fPIC -fno-range-check -O
#  -fcheck=all

FC=gfortran

FFLAGS += -I../lib
LDFLAGS = -L../lib
LDLIBS  = -lhe3

all: libtext1r.a libtext1r.so octave doc

doc:
	make -C doc

libtext1r.a: text1r.o ../lib/libhe3.a ../lib/external/libtn.a
	ar rs $@ $+

libtext1r.so: text1r.o ../lib/libhe3.a ../lib/external/libtn.a
	$(FC) --shared -fPIC -o $@ $+

text1r.fh text1r.h: text1r.def make_inc
	./make_inc

text1r.o: text1r.f90 text1r.fh
	$(FC) -c $< $(FFLAGS)

octave: libtext1r.so
	rm -f *.mex
	octave -q --eval 'mex_text1r'

clean:
	rm -f *.o *.a *.so
	rm -f text1r.fh text1r.h
	rm -f *.mexglx *.mex
