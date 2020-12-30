#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

//random float between -1 and +1
float RandomNumber()
{
  return (float)rand()/(float)RAND_MAX * 2.0f - 1.0;
}

int main(int argc, char *argv[])
{
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  srand(time(NULL)+rank);
  // 'parsing' optional input parameter = problem size
  long SIZE = 20000;
  if (argc > 1) {
    SIZE = atol(argv[1]);
  }
  long in = 0;
  long out = 0;
  long reducedin = 0;
  long reducedout = 0;
  long reducedsize = 0;
  for (long i = 0; i < SIZE; ++i) {
    float x = RandomNumber();
    float y = RandomNumber();
    ((x * x + y * y) < 1.0f) ? in++ : out++;
  }

  printf("inside=%ld outside=%ld total=%ld \n", in ,out, SIZE);
  printf("%f\n", (double)in/(double)SIZE*4.0);

  MPI_Reduce(&in, &reducedin, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&out, &reducedout, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&SIZE, &reducedsize, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0){
    printf("final inside=%ld outside=%ld total=%ld \n", reducedin ,reducedout, reducedsize);
    printf("%f\n", (double)reducedin/(double)reducedsize*4.0);
  }
  
  MPI_Finalize();
  return EXIT_SUCCESS;
}
