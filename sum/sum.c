#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <linux/time.h>
#include <mm_malloc.h>
#include "../nbody/Common/timers.h"

#define NITER 10000


void clobber() {
  __asm__ __volatile__ ("" : : : "memory");
}

static void escape(void* p)
{
    __asm__ __volatile__("" : : "g"(p) : "memory");
}

#define N 1000000

float sum(const float* __restrict__ array){
	float sum = 0.f;
#pragma omp parallel for shared(array) reduction(+: sum)
//#pragma unroll 
	for (int i = 0; i < N; i += 1) {
		sum += array[i+0];
#if 0
		sum += array[i+1];
		sum += array[i+2];
		sum += array[i+3];
		sum += array[i+4];
		sum += array[i+5];
		sum += array[i+6];
		sum += array[i+7];
		sum += array[i+8];
		sum += array[i+9];
		sum += array[i+10];
		sum += array[i+11];
		sum += array[i+12];
		sum += array[i+12];
		sum += array[i+14];
#endif
	}
	return sum;
}
//
float sum1(const float* __restrict__ array){
        float sum1 = 0.f;
        float sum2 = 0.f;
#pragma omp parallel for shared(sum1, sum2, array) reduction(+: sum1, sum2)
        for (int i = 0; i < N; i += 8) {
                sum1 += array[i+0];
                sum2 += array[i+1];
#if 1
                sum1 += array[i+2];
                sum2 += array[i+3];
                sum1 += array[i+4];
                sum2 += array[i+5];
                sum1 += array[i+6];
                sum2 += array[i+7];
#endif
        }
        return sum1 + sum2;
}
//
float sum2(const float* __restrict__ array){
        float sum1 = 0.f;
        float sum2 = 0.f;
        float sum3 = 0.f;
        float sum4 = 0.f;
#pragma omp parallel for shared(array) reduction(+: sum1, sum2, sum3, sum4)
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
#pragma omp parallel for shared(array) reduction(+: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8)
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
        return ((sum1 + sum2) + (sum3 + sum4) + (sum5 + sum6) + (sum7 + sum8));
}
//
float sum4(const float* __restrict__ array){
        float sum1  = 0.f;
        float sum2  = 0.f;
        float sum3  = 0.f;
        float sum4  = 0.f;
        float sum5  = 0.f;
        float sum6  = 0.f;
        float sum7  = 0.f;
        float sum8  = 0.f;
        float sum9  = 0.f;
        float sum10 = 0.f;
        float sum11 = 0.f;
        float sum12 = 0.f;
        float sum13 = 0.f;
        float sum14 = 0.f;
        float sum15 = 0.f;
        float sum16 = 0.f;
#pragma omp parallel for shared(array) reduction(+: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16)
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
                sum8 += array[i+8];
        }
        return ((sum1 + sum2) + (sum3 + sum4)) + ((sum5 + sum6) + (sum7 + sum8)) + sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16;
}
//
//
int main(const int argc, const char** argv)
{
	printf("Buffer size = %f MB\n", N/1024/1024.);
	//
	float* array = (float*) _mm_malloc(N*sizeof(float), 64);
	for (int ii = 0; ii < NITER; ii += 1)
		for (int jj = 0; jj < N; jj += 1) array[jj] = (float) rand()/(float) RAND_MAX;
	//
	//
	double time;
	//
#if 0
	{
		float mysum = 0.;
		time = -myseconds(); 
		unsigned long c0 = rdpmc_actual_cycles();
		for (int ii = 0; ii < NITER; ii += 1)
		{
			clobber();
			mysum += sum(array);
		}
		unsigned long c1 = rdpmc_actual_cycles();
		time += myseconds();
		printf("sum0 = %.15f, throughput = %f GB/s [%f s], %f GF/s\n", (double) mysum, sizeof(float)*N*(NITER/time/1024/1024/1024), time,  N*NITER/time/1e9);
	}
#endif
	//
#if 0
	{
		float mysum = 0.;
		time = -myseconds();
		unsigned long c0 = rdpmc_actual_cycles();
		for (int ii = 0; ii < NITER; ii += 1)
		{
			clobber();
			mysum += sum1(array);
		}
		unsigned long c1 = rdpmc_actual_cycles();
		time += myseconds();
		printf("sum1 = %.15f, throughput = %f GB/s [%f s], %f GF/s\n", (double) mysum, sizeof(float)*N*(NITER/time/1024/1024/1024), time, N*NITER/time/1e9);
	}
#endif
	//
#if 0
	{
                float mysum = 0.;
                time = -myseconds();
		unsigned long c0 = rdpmc_actual_cycles();
                for (int ii = 0; ii < NITER; ii += 1)
		{
			clobber();
                        mysum += sum2(array);
		}
		unsigned long c1 = rdpmc_actual_cycles();
                time += myseconds();
                printf("sum2 = %.15f, throughput = %f GB/s [%f s], %f GF/s\n", (double) mysum, sizeof(float)*N*(NITER/time/1024/1024/1024), time, N*NITER/time/1e9);
        }	
#endif
	//
#if 0
	{
                float mysum = 0.;
                time = -myseconds();
                for (int ii = 0; ii < NITER; ii += 1)
		{
			clobber();
                        mysum += sum3(array);
		}
                time += myseconds();
                printf("sum3 = %.15f, throughput = %f GB/s [%f s], %f GF/s\n", (double) mysum, sizeof(float)*N*(NITER/time/1024/1024/1024), time, N*(NITER/time)/1.e9);
        }	
#endif
	//
	_mm_free(array);
}
