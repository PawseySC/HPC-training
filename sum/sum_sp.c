#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <linux/time.h>
#include <mm_malloc.h>

#define NITER 10000

double myseconds()
{
        struct timespec mysec;
        clock_gettime( CLOCK_REALTIME, &mysec);
        return ( (double) mysec.tv_sec + (double) mysec.tv_nsec * 1.e-9 );
}


void clobber() {
  __asm__ __volatile__ ("" : : : "memory");
}

static void escape(void* p)
{
    __asm__ __volatile__("" : : "g"(p) : "memory");
}

#define N 100000

float sum(const float* __restrict__ array){
	float sum = 0.f;
#pragma omp parallel for shared(sum, array) reduction(+: sum)
	for (int i = 0; i < N; i += 8) {
		sum += array[i+0];
		sum += array[i+1];
		sum += array[i+2];
		sum += array[i+3];
		sum += array[i+4];
		sum += array[i+5];
		sum += array[i+6];
		sum += array[i+7];
	}
	return sum;
}
//
float sum1(const float* __restrict__ array){
        float sum1 = 0.f;
        float sum2 = 0.f;
#pragma omp parallel for shared(sum, array) reduction(+: sum1, sum2)
        for (int i = 0; i < N; i += 8) {
                sum1 += array[i+0];
                sum2 += array[i+1];
                sum1 += array[i+2];
                sum2 += array[i+3];
                sum1 += array[i+4];
                sum2 += array[i+5];
                sum1 += array[i+6];
                sum2 += array[i+7];
        }
        return sum1 + sum2;
}
//
float sum2(const float* __restrict__ array){
        float sum1 = 0.f;
        float sum2 = 0.f;
        float sum3 = 0.f;
        float sum4 = 0.f;
#pragma omp parallel for shared(sum, array) reduction(+: sum1, sum2, sum3, sum4)
        for (int i = 0; i < N; i += 8) 
	{
                sum1 += array[i+0];
                sum2 += array[i+1];
                sum3 += array[i+2];
                sum4 += array[i+3];
                sum1 += array[i+4];
                sum2 += array[i+5];
                sum3 += array[i+6];
                sum4 += array[i+7];
        }
        return (sum1 + sum2) + (sum3 + sum4);
}
//
float sum3(const float* __restrict__ array){
        float sum1 = 0.f;
        float sum2 = 0.f;
        float sum3 = 0.f;
        float sum4 = 0.f;
        float sum5 = 0.f;
        float sum6 = 0.f;
        float sum7 = 0.f;
        float sum8 = 0.f;
#pragma omp parallel for shared(sum, array) reduction(+: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8)
        for (int i = 0; i < N; i += 8) 
	{
                sum1 += array[i+0];
                sum2 += array[i+1];
                sum3 += array[i+2];
                sum4 += array[i+3];
                sum5 += array[i+4];
                sum6 += array[i+5];
                sum7 += array[i+6];
                sum8 += array[i+7];
        }
        return ((sum1 + sum2) + (sum3 + sum4)) + ((sum5 + sum6) + (sum7 + sum8));
}
//
//
int main(const int argc, const char** argv)
{
	printf("Buffer size = %f MB\n", N/1024/1024.);
	//
	float* array = (float*) _mm_malloc(N*sizeof(float), 64);
	for (int ii = 0; ii < N; ii += 1) array[ii] = (float) rand()/(float) RAND_MAX;
	//
	//
	double time;
	//
	{
		float mysum = 0.;
		time = -myseconds(); 
		for (int ii = 0; ii < NITER; ii += 1)
		{
			clobber();
			mysum += sum(array);
		}
		time += myseconds();
		printf("sum0 = %.15f, throughput = %f GB/s [%f s]\n", (double) mysum, N*(NITER/time/1024/1024/1024), time);
	}
	//
	{
		float mysum = 0.;
		time = -myseconds();
		for (int ii = 0; ii < NITER; ii += 1)
		{
			clobber();
			mysum += sum1(array);
		}
		time += myseconds();
		printf("sum1 = %.15f, throughput = %f GB/s [%f s]\n", (double) mysum, N*(NITER/time/1024/1024/1024), time);
	}
	//
	{
                float mysum = 0.;
                time = -myseconds();
                for (int ii = 0; ii < NITER; ii += 1)
		{
			clobber();
                        mysum += sum2(array);
		}
                time += myseconds();
                printf("sum2 = %.15f, throughput = %f GB/s [%f s]\n", (double) mysum, N*(NITER/time/1024/1024/1024), time);
        }	
	//
	{
                float mysum = 0.;
                time = -myseconds();
                for (int ii = 0; ii < NITER; ii += 1)
		{
			clobber();
                        mysum += sum3(array);
		}
                time += myseconds();
                printf("sum3 = %.15f, throughput = %f GB/s [%f s]\n", (double) mysum, N*(NITER/time/1024/1024/1024), time);
        }	
	//
	_mm_free(array);
}
