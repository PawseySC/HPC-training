#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <linux/time.h>

static
__inline__ uint64_t rdtsc(void) {
        uint32_t l, h;
        __asm__ __volatile__ ("xorl %%eax,%%eax\ncpuid"
                              ::: "%rax", "%rbx", "%rcx", "%rdx");
        __asm__ __volatile__ ("rdtsc" : "=a" (l), "=d" (h));
        return (uint64_t)h << 32 | l;
}

static
__inline__
unsigned long long rdtscp()
{
        unsigned int cycles_high, cycles_low;
        __asm__ volatile (
                        "RDTSCP\n\t"/*read the clock*/
                        "mov %%edx, %0\n\t"
                        "mov %%eax, %1\n\t"
                        "CPUID\n\t": "=r" (cycles_high), "=r"
                        (cycles_low):: "%rax", "%rbx", "%rcx", "%rdx");
        return ((unsigned long long)cycles_low + ((unsigned long long)cycles_high << 32));
}

unsigned long rdpmc_instructions()
{
   unsigned a, d, c;

   c = (1<<30);
   __asm__ volatile("rdpmc" : "=a" (a), "=d" (d) : "c" (c));

   return ((unsigned long)a) | (((unsigned long)d) << 32);;
}

unsigned long rdpmc_actual_cycles()
{
   unsigned a, d, c;

   c = (1<<30)+1;
   __asm__ volatile("rdpmc" : "=a" (a), "=d" (d) : "c" (c));

   return ((unsigned long)a) | (((unsigned long)d) << 32);;
}

unsigned long rdpmc_reference_cycles()
{
   unsigned a, d, c;

   c = (1<<30)+2;
   __asm__ volatile("rdpmc" : "=a" (a), "=d" (d) : "c" (c));

   return ((unsigned long)a) | (((unsigned long)d) << 32);;
}

// needs -std=c11 to work
double myseconds()
{
        struct timespec mysec;
        clock_gettime( CLOCK_REALTIME, &mysec);
        return ( (double) mysec.tv_sec + (double) mysec.tv_nsec * 1.e-9 );
}

