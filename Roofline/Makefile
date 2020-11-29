CC = gcc
CFLAGS = -O3 -march=skylake-avx512 -DSTREAM_ARRAY_SIZE=80000000 -DNTIMES=20 
CFLAGS += -fopenmp 
CFLAGS += -lopenblas

FF = g77
FFLAGS = -O3

all: stream dgemm

stream: stream.c
	$(CC) $(CFLAGS) stream.c -o stream

dgemm: dgemm.c
	$(CC) $(CFLAGS) dgemm.c -o dgemm 

clean:
	rm -f stream dgemm

# an example of a more complex build line for the Intel icc compiler
stream.icc: stream.c
	cc -ffreestanding -qopenmp -DSTREAM_ARRAY_SIZE=80000000 -DNTIMES=20 stream.c -o stream.omp.AVX2.80M.20x.icc
