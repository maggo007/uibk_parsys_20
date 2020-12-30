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
  int SIZE = 20000;
  if (argc > 1) {
    SIZE = atoi(argv[1]);
  }
  int in = 0;
  int out = 0;
  int numthreads;
  numthreads = omp_get_num_threads();
  int chunk = SIZE/numthreads;
  float x,y;
#pragma omp parallel for private(x,y) reduction(+:in,out)
  for (int i = 0; i < SIZE; ++i) {
    unsigned int seed = omp_get_thread_num();
    x = RandomNumber(&seed);
    y = RandomNumber(&seed);
    ((x * x + y * y) < 1.0f) ? in++ : out++;
  }
  

  printf("inside=%d outside=%d total=%d numberofthreads=%d\n", in ,out, SIZE, numthreads);
  printf("%f\n", (double)in/(double)(SIZE)*4.0);
  if (in + out != SIZE){
    printf("sum of in and out not adding up to SIZE\n");
  }
  return EXIT_SUCCESS;
}
