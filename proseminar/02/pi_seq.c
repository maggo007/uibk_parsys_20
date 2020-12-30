#include <stdio.h>
#include <stdlib.h>

//random float between -1 and +1
float RandomNumber()
{
  return (float)rand()/(float)RAND_MAX * 2.0f - 1.0;
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
  for (long i = 0; i < SIZE; ++i) {
    float x = RandomNumber();
    float y = RandomNumber();
    ((x * x + y * y) < 1.0f) ? in++ : out++;
  }

  printf("inside=%ld outside=%ld total=%ld \n", in ,out, SIZE);
  printf("%f\n", (double)in/(double)SIZE*4.0);

  return EXIT_SUCCESS;
}
