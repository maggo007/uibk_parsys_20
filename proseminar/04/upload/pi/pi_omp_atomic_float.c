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
  float pi = 0;
  long numthreads;
#pragma omp parallel
  {
    long localin = 0;
    long localout = 0;
    float localpi = 0;
    numthreads = omp_get_num_threads();
    long chunk = SIZE/numthreads;
    srand(time(NULL)+omp_get_thread_num());
    unsigned int seed = omp_get_thread_num();
    for (long i = 0; i < chunk; ++i) {
      float x = RandomNumber(&seed);
      float y = RandomNumber(&seed);
      ((x * x + y * y) < 1.0f) ? localin++ : localout++;
    }
    localpi = (float)localin/(float)chunk*4.0;
#pragma omp atomic
    pi +=localpi/numthreads;
  }

  printf("numberofthreads=%ld\n", numthreads);
  printf("%f\n", pi);
  return EXIT_SUCCESS;
}
