#pragma once

#define NUM_BODIES  32768
#define DT          0.01
#define NITERS      10


#if !defined(__SP__) && !defined(__DP__)
#define __SP__
#endif

#ifdef __SP__
#define TILE_SIZE   16
#warning "single precission"
#define SOFTENING   1e-20
#define FLOAT float
#define MPIFLOAT MPI_FLOAT
#define INVSQRT  1.f/sqrtf
#define SQRT  sqrtf
#define ZERO  0.f
#define ONE   1.f
#define TWO   2.f
#define THREE 3.f
#endif

#ifdef __DP__ 

#define TILE_SIZE 8
inline
double marla_inv_sqrt(double ix) 	
{
	double w = ix*0.5;
	double x = ix;
	x = (double) (1.f/sqrtf((float) x));
        x = x*(1.5 - w*x*x);
	//x = x*(1.5 - w*x*x);
	return x;
}


#warning "double precision"
#define SOFTENING   1e-20
#define FLOAT double
#define INVSQRT marla_inv_sqrt
#define SQRT sqrt
#define ZERO  0.
#define ONE   1.
#define TWO   2.
#define THREE 3.
#endif

