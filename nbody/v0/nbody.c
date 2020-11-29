#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "timers.h"
#include "defs.h"

typedef struct 
{ 
	FLOAT x; 
	FLOAT y; 
	FLOAT z; 
	FLOAT vx; 
	FLOAT vy; 
	FLOAT vz; 
} Body;


void setUp(Body* __restrict__ body, const int n)
{
        for (int i = 0; i < n; i++)
        {
                body[i].x  = TWO*((FLOAT) rand()/(FLOAT)RAND_MAX) - ONE;
                body[i].y  = TWO*((FLOAT) rand()/(FLOAT)RAND_MAX) - ONE;
                body[i].z  = TWO*((FLOAT) rand()/(FLOAT)RAND_MAX) - ONE;
                body[i].vx = TWO*((FLOAT) rand()/(FLOAT)RAND_MAX) - ONE;
                body[i].vy = TWO*((FLOAT) rand()/(FLOAT)RAND_MAX) - ONE;
                body[i].vz = TWO*((FLOAT) rand()/(FLOAT)RAND_MAX) - ONE;
        }
}


void updatePos(Body* __restrict__ body, const FLOAT dt, const int n) 
{
	for (int i = 0 ; i < n; i++)
	{ // integrate position
		body[i].x += body[i].vx*dt;
		body[i].y += body[i].vy*dt;
		body[i].z += body[i].vz*dt;
	}
}


FLOAT energy(Body *p, const int n)
{
	FLOAT enr = ZERO;
	for (int i = 0; i < n; i++)
	{
		enr += p[i].vx*p[i].vx;
		enr += p[i].vy*p[i].vy;
		enr += p[i].vz*p[i].vz;
	}
	return enr/2.;
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
		for (int j = 0; j < n; j++) 
		{
			FLOAT dx = p[j].x - p[i].x; // 1 flop 2 loads
			FLOAT dy = p[j].y - p[i].y; // 1 flop 2 loads
			FLOAT dz = p[j].z - p[i].z; // 1 flop 2 loads
			//
			FLOAT dist        = dx*dx + dy*dy + dz*dz + SOFTENING; // 6 flops
			FLOAT distPow32 = pow(dist, 3./2.);                    // 1 fLops? 
			//
			Fx += dx/distPow32; // 2 fLops
			Fy += dy/distPow32; // 2 flops
			Fz += dz/distPow32; // 2 flops
		}  // ~16 flops and 6 loads per iteration  
		//
		p[i].vx += dt*Fx; // 2 flops
		p[i].vy += dt*Fy; // 2 flops
		p[i].vz += dt*Fz; // 2 flops
	}
}

int main(const int argc, const char** argv) 
{
	int nBodies = NUM_BODIES;
	if (argc > 1) nBodies = atoi(argv[1]);
	//
	const FLOAT dt  = 0.01; // time step
	const int nIters = 10;  // simulation iterations
	//
	int body_size = nBodies*sizeof(Body);
	void *buf   = (void*) malloc(body_size);
	Body *p       = (Body*)buf;
	//
	setUp(p, nBodies); // Init pos / vel data
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
			bodyForce(p, dt, nBodies); // compute interbody forces
			updatePos(p, dt, nBodies);
		unsigned long c1 = rdpmc_actual_cycles();
		time += myseconds();
		totalTime += time;
		//
		printf("%4d %10.3e %10.3e    \033[42m%.1f\033[0m    %3.1f\n", iter, time,  (double) nBodies*nBodies/time, 21.*nBodies*nBodies*1e-9/time, 30.*nBodies*nBodies/((double) c1 - c0));
	}
	double avgTime = totalTime / (double)(nIters-1); 

        FLOAT en = energy(p, nBodies);
        printf("\n Average performance = %.5f, energy = %g\n",  21.*nBodies*nBodies*1e-9/avgTime, en);
}
