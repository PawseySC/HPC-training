CC=gcc

CFLAGS += -std=c11 -g -Ofast
CFLAGS += -I../Common

ifeq ($(CC), gcc)
	CFLAGS += -march=skylake-avx512
	CFLAGS += -mprefer-vector-width=512
	CFLAGS += -fopenmp
	CFLAGS += -lm
else
	CFLAGS += -xCORE-AVX512
	CFLAGS += -qopt-zmm-usage=high
	CFLAGS += -qopenmp
endif



ADVPRJ = "./adv"

all:
	#$(CC) nbody.c -o nbody $(CFLAGS) 
	$(CC) nbody.c -o nbody_sp $(CFLAGS) -D__SP__
	#
	#$(CC) nbody.c -o nbody $(CFLAGS) 
	$(CC) nbody.c -o nbody_dp $(CFLAGS) -D__DP__
	# 

survey:
	advixe-cl -collect survey -project-dir $(ADVPRJ) -- ./nbody_sp 65536
	advixe-cl -c tripcounts -flop -project-dir $(ADVPRJ) --search-dir src:=./ --search-dir bin:=./ -- ./nbody_sp 65536
	

roofline:
	#advixe-cl -collect survey -project-dir $(ADVPRJ) -- ./nbody_sp 65536
	#advixe-cl -collect tripcounts -flop -project-dir $(ADVPRJ) -- ./nbody_sp 65536
	#advixe-cl -collect roofline -project-dir $(ADVPRJ) --search-dir src:=./ --search-dir bin:=./ -- ./nbody_sp 65536
	advixe-cl --collect roofline --integrated -project-dir $(ADVPRJ) --search-dir src:=./ --search-dir bin:=./ -- ./nbody_sp 65536
	advixe-cl --snapshot -project-dir $(ADVPRJ) --pack --cache-sources --cache-binaries --search-dir src:=./ --search-dir bin:=./ -- ./nbody_sp 65536

clean:
	rm *.optrpt nbody_sp nbody_dp
