all: matvec

CXXFLAGS = -Wall -g -O2 -march=native -fopenmp
LDFLAGS = -lflexiblas -fopenmp

matvec: matvec_main.o matvec_func.o
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	rm -f matvec_main.o matvec_func.o matvec
