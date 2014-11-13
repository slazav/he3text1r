FFLAGS= -Werror -Wconversion\
  -Wline-truncation\
  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\
  -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow\
  -Wno-unused-parameter -fPIC -fno-range-check -O
#  -fcheck=all

FC=gfortran

HE3LIB_PATH = ../he3lib

FFLAGS += -I$(HE3LIB_PATH)
LDFLAGS = -L$(HE3LIB_PATH)
LDLIBS  = -lhe3

all: libtext1r.a libtext1r.so octave doc


doc:
	make -C doc

text1r.fh text1r.h: text1r.def make_text1r
	./make_text1r

libtext1r.a: text1r.o tn/tn.a
	ar rs $@ $+

libtext1r.so: text1r.o tn/tn.a
	$(FC) --shared -fPIC -o $@ $+

tn/tn.a:
	make -C tn

text1r.o: text1r.f90 text1r.fh
	$(FC) -c $< $(FFLAGS)

octave: libtext1r.so libhe3.so
	rm -f *.mex
	octave -q --eval 'mex_text1r'

libhe3.so: $(HE3LIB_PATH)/libhe3.so
	ln -s $(HE3LIB_PATH)/libhe3.so .

clean:
	rm -f *.o *.a *.so
	rm -f text1r.fh text1r.h
	rm -f *.mexglx *.mex
