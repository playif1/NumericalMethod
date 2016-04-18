#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <signal.h>
#include <fenv.h>
#include <math.h>


void handler(int sig, siginfo_t *sip)
{
   // print out the exception type and it's address
   fprintf(stderr, "fp exception %x at address %x \n", sig, (unsigned) sip);
   exit(1);
}

int main() {
   double x;
   /* trap on common floating point exceptions */
   // set the floating pointing exception signal for handler
   void (*prev_handler)(int);
   prev_handler = signal(SIGFPE, handler);

   // enable the "common" floating point exception
   feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

   /* cause an underflow exception (not reported) */
   x = DBL_MIN;
   printf("min_normal = %g \n", x);
   x = x / 13.0;
   printf("min_normal / 13.0 = %g\n", x);

   /* cause an overflow exception (reported) */
   x = DBL_MAX;
   printf("max_normal = %g\n", x);
   x = x * x;
   printf("max_normal * max_normal = %g\n",x);
   return 0;
}
