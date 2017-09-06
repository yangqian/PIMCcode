LDFLAGS= -openmp -lgsl -lgslcblas -lm
CC=mpic++
FERMIFLAGS=-Wall -O2 #-Wliteral-suffix# -Wno-write-strings -Dll #ll flag for long double
SOURCES=$(wildcard ./*.c)

OBJECTS=$(patsubst ./%,fermi/%,$(patsubst %.c,%.o,$(wildcard ./*.c)))
EXECUTABLE=run

all: $(SOURCES) ./fermi/$(EXECUTABLE)

./fermi/$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

./fermi/%.o:%.c 
	$(CC) -c $(FERMIFLAGS) $< -o $@

$(OBJECTS) : ext.h func.h myran.h
./fermi/main.o : var.h
./fermi/DoubGridEstimator.o : LinearGridEstimator.h accumulator.h
./fermi/LinearGridEstimator.o : LinearGridEstimator.h accumulator.h
./fermi/calc.o : LinearGridEstimator.h accumulator.h dilation.h moves.h
./fermi/dilation.o : dilation.h zerorange.h
./fermi/initialize.o : moves.h

.PHONY: clean
clean:
	rm -rf ./fermi/*o ./fermi/test ./test/*.dat ./test/3/*.dat ./test/fpath.save ./test/3/fpath.save ./fermi/run
