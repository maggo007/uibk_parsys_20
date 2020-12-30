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
#pragma omp parallel for reduction(+:in,out)
  for (long i = 0; i < SIZE; ++i) {
    unsigned int seed = omp_get_thread_num()+i*i;
    float x = RandomNumber(&seed);
    float y = RandomNumber(&seed);
    ((x * x + y * y) < 1.0f) ? in++ : out++;
  }

  printf("inside=%ld outside=%ld total=%ld\n", in ,out, SIZE);
  printf("%f\n", (double)in/(double)(SIZE)*4.0);
  if (in + out != SIZE){
    printf("sum of in and out not adding up to SIZE\n");
  }
  return EXIT_SUCCESS;
}
