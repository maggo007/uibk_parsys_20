#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

//random float between -1 and +1
float RandomNumber(unsigned int* seed)
{
  return (float)rand_r(seed)/(float)RAND_MAX * 2.0f - 1.0;
}

int main(int argc, char *argv[])
{
  // 'parsing' optional input parameter = problem size
  long SIZE = 20000;
  if (argc > 1) {
    SIZE = atol(argv[1]);
  }
  long in = 0;
  long out = 0;
  long numthreads;
#pragma omp parallel reduction(+:in,out)
  {
    long localin = 0;
    long localout = 0;
    numthreads = omp_get_num_threads();
    long chunk = SIZE/numthreads;
    srand(time(NULL)+omp_get_thread_num());
    unsigned int seed = omp_get_thread_num();
    float x,y;
    for (long i = 0; i < chunk; ++i) {
      x = RandomNumber(&seed);
      y = RandomNumber(&seed);
      ((x * x + y * y) < 1.0f) ? localin++ : localout++;
    }
    in +=localin;
    out +=localout;
  }

  printf("inside=%ld outside=%ld total=%ld numberofthreads=%ld\n", in ,out, SIZE, numthreads);
  printf("%f\n", (double)in/(double)(SIZE)*4.0);
  if (in + out != SIZE){
    printf("sum of in and out not adding up to SIZE\n");
  }
  return EXIT_SUCCESS;
}
