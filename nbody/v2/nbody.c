#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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
#pragma omp parallel for schedule(guided) 
        for (int i = 0; i < n; i++)
        {
                FLOAT Fx = 0.;
                FLOAT Fy = 0.;
                FLOAT Fz = 0.;
                //
                #pragma simd
		#pragma unroll
                for (int j = 0; j < n; j++)
                {
                        const FLOAT dx = p->x[j] - p->x[i]; // 1 flop 2 loads
                        const FLOAT dy = p->y[j] - p->y[i]; // 1 flop 2 loads
                        const FLOAT dz = p->z[j] - p->z[i]; // 1 flop 2 loads
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
                p->vx[i] += dt*Fx; // 2 flops
                p->vy[i] += dt*Fy; // 2 flops
                p->vz[i] += dt*Fz; // 2 flops
         }
}


int main(const int argc, const char** argv) 
{
	int nBodies = NUM_BODIES;
	if (argc > 1) nBodies = atoi(argv[1]);

	const FLOAT dt   = DT; // time step
	const int nIters = NITERS;  // simulation iterations

	int body_size = nBodies*sizeof(FLOAT);
	//
	printf("\n");
	printf("\n\033[1mNBody: num particles = %d, dt = %.3f, sizeof(FLOAT) = %d\033[0m\n", nBodies, dt, sizeof(FLOAT));
	printf("\n");
	//
	Body p;       
	//
	p.x  = (FLOAT*) malloc(body_size); 
	p.y  = (FLOAT*) malloc(body_size); 
	p.z  = (FLOAT*) malloc(body_size); 
	//
	p.vx = (FLOAT*) malloc(body_size); 
	p.vy = (FLOAT*) malloc(body_size); 
	p.vz = (FLOAT*) malloc(body_size); 
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
		}
		unsigned long c1 = rdpmc_actual_cycles();
		updatePos(&p, dt, nBodies);
		//
		time += myseconds();
		totalTime += time;
		//
		printf("%4d %10.3e %10.3e    \033[42m%.1f\033[0m    %3.1f\n", iter, time,  (double) nBodies*nBodies/time, 21.*nBodies*nBodies*1e-9/time, 30.*nBodies*nBodies/((double) c1 - c0));
	}
	double avgTime = totalTime / (FLOAT)nIters;

	FLOAT en = energy(&p, nBodies);
	 printf("\n\033[1m Average performance = \033[42m%.1f\033[0m, energy = %g\033[0m\n\n",  21.*nBodies*nBodies*1e-9/avgTime, en);

	free(p.x);
	free(p.y);
	free(p.z);
	//
	free(p.vx);
	free(p.vy);
	free(p.vz);
}
