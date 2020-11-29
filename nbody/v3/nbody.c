#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mm_malloc.h>
#include "timers.h"
#include "defs.h"

typedef struct 
{ 
	FLOAT* __restrict__ x; 
	FLOAT* __restrict__ y; 
	FLOAT* __restrict__ z; 
	FLOAT* __restrict__ vx; 
	FLOAT* __restrict__ vy; 
	FLOAT* __restrict__ vz; 
} Body;


void setUp(Body* __restrict__ body, const int n) 
{
	printf("%f\n", (FLOAT) RAND_MAX);
	for (int i = 0; i < n; i++) 
	{
		body->x[i]  = TWO*((FLOAT) rand()/(FLOAT)RAND_MAX) - ONE;
		body->y[i]  = TWO*((FLOAT) rand()/(FLOAT)RAND_MAX) - ONE;
		body->z[i]  = TWO*((FLOAT) rand()/(FLOAT)RAND_MAX) - ONE;
		body->vx[i] = TWO*((FLOAT) rand()/(FLOAT)RAND_MAX) - ONE;
		body->vy[i] = TWO*((FLOAT) rand()/(FLOAT)RAND_MAX) - ONE;
		body->vz[i] = TWO*((FLOAT) rand()/(FLOAT)RAND_MAX) - ONE;
	}
}


void updatePos(Body* __restrict__ body, const FLOAT dt, const int n) 
{
#pragma omp parallel for schedule(guided)
	for (int i = 0 ; i < n; i++)
	{ // integrate position
		body->x[i] += body->vx[i]*dt;
		body->y[i] += body->vy[i]*dt;
		body->z[i] += body->vz[i]*dt;
	}
}


FLOAT energy(Body* __restrict__ body, const int n)
{
	FLOAT enr = ZERO;
	for (int i = 0; i < n; i++)
	{
		enr += body->vx[i]*body->vx[i];
		enr += body->vy[i]*body->vy[i];
		enr += body->vz[i]*body->vz[i];
	}
	return enr/TWO;
}




// the main routine
void bodyForce(Body *p, FLOAT const dt, int n) 
{
	//long int count = 0;
#pragma omp parallel for schedule(guided) 
	for (int ii = 0; ii < n; ii += 1) 
	{ 
		FLOAT Fx = ZERO;
		FLOAT Fy = ZERO;
		FLOAT Fz = ZERO;
		//
		#pragma simd
		#pragma unroll
                for (int jj = 0; jj < n; jj++)
                {
                        const FLOAT dx = p->x[jj] - p->x[ii]; // 1 flop 2 loads
                        const FLOAT dy = p->y[jj] - p->y[ii]; // 1 flop 2 loads
                        const FLOAT dz = p->z[jj] - p->z[ii]; // 1 flop 2 loads
                        //
                        const FLOAT dist        = dx*dx + dy*dy + dz*dz + (FLOAT) SOFTENING; // 6 flops
                        const FLOAT invSqrtDist = ONE/SQRT(dist);                     // 2 fLops (?)
                        //
                        const FLOAT invDist3 = invSqrtDist*invSqrtDist*invSqrtDist;  // 3 flops
                        //
                        Fx += dx*invDist3; // 2 fLops
                        Fy += dy*invDist3; // 2 flops
                        Fz += dz*invDist3; // 2 flops
                }  // 20 flops,
                //
                p->vx[ii] += dt*Fx; // 2 flops
                p->vy[ii] += dt*Fy; // 2 flops
                p->vz[ii] += dt*Fz; // 2 flops
	}
}



int main(const int argc, const char** argv) 
{
	int nBodies = NUM_BODIES;
	if (argc > 1) nBodies = atoi(argv[1]);
	//
	const FLOAT dt   = DT; // time step
	const int nIters = NITERS;  // simulation iterations
	//
	int body_size = nBodies*sizeof(FLOAT);
	//
	Body p;       //= (Body*) buf;
	//
	p.x  = (FLOAT*) _mm_malloc(body_size, 64); 
	p.y  = (FLOAT*) _mm_malloc(body_size, 64); 
	p.z  = (FLOAT*) _mm_malloc(body_size, 64); 
	//
	p.vx = (FLOAT*) _mm_malloc(body_size, 64); 
	p.vy = (FLOAT*) _mm_malloc(body_size, 64); 
	p.vz = (FLOAT*) _mm_malloc(body_size, 64); 
	//
	setUp(&p, nBodies); // Init pos / vel data
	//
        printf("Iter   Time [s]    Inter/s  Gflops flops/c\n");
        printf("------------------------------------------\n");
	//
	double totalTime = 0.;
	//
	for (int iter = 1; iter <= nIters; iter++) 
	{
		double time = -myseconds();
		unsigned long c0 = rdpmc_actual_cycles();
		{
			bodyForce(&p, dt, nBodies); // compute interbody forces
			updatePos(&p, dt, nBodies);
		}
		unsigned long c1 = rdpmc_actual_cycles();
		time += myseconds();
		totalTime += time;
		printf("%4d %10.3e %10.3e    \033[42m%.1f\033[0m    %3.1f\n", iter, time,  (double) nBodies*nBodies/time, 21.*nBodies*nBodies*1e-9/time, 30.*nBodies*nBodies/((double) c1 - c0));
	}
	double avgTime = totalTime / (FLOAT)nIters;

	FLOAT en = energy(&p, nBodies);
	//printf("\n Average performance = %.5f, energy = %g\n",  21.*nBodies*nBodies*1e-9/avgTime, en);
	printf("\n\033[1m Average performance = \033[42m%.1f\033[0m, energy = %g\033[0m\n\n",  21.*nBodies*nBodies*1e-9/avgTime, en);
	//
	_mm_free(p.x);
	_mm_free(p.y);
	_mm_free(p.z);
	//
	_mm_free(p.vx);
	_mm_free(p.vy);
	_mm_free(p.vz);
	//
}
