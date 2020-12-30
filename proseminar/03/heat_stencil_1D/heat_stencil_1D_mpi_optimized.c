#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

typedef double value_t;

#define RESOLUTION 120

// -- vector utilities --

typedef value_t *Vector;

Vector createVector(int N);

void releaseVector(Vector m);

void printTemperature(Vector m, int N);

// -- simulation code ---

int main(int argc, char **argv) {
  int rank;
  int size;
  int partsize;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  printf("rank=%d with size=%d\n", rank, size);
  // 'parsing' optional input parameter = problem size
  int N = 2000;
  if (argc > 1) {
    N = atoi(argv[1]);
  }
  int T = N * 500;
  printf("Computing heat-distribution for room size N=%d for T=%d timesteps\n",
         N, T);
  if (N % size != 0) {
    printf("size not ok for problem size. not dividable without remainder");
    exit(EXIT_FAILURE);
  }
  partsize = N / size;

  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Vector A = createVector(N);

  // set up initial conditions in A
  for (int i = 0; i < N; i++) {
    A[i] = 273; // temperature is 0° C everywhere (273 K)
  }

  // and there is a heat source in one corner
  int source_x = N / 4;
  A[source_x] = 273 + 60;

  printf("Initial:\t");
  printTemperature(A, N);
  printf("\n");

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(N);

  // for each time step ..
  long long start = rank * partsize;
  long long stop = rank * partsize + partsize;
  for (int t = 0; t < T; t++) {
    // make two waves to send and receive overlapping parts (now only using one)
    // .. we propagate the temperature for the assigned part of N[(rank*partsize)+i,(rank*partsize)+(i+1),...,(ranke*partsize)+(partsize-1)] for i < partsize
    // example for N = 12, partsize=4, size=3, rank=0, operates on N[(0*4)+0, (0*4)+1, (0*4)+2, (0*4)*3] --> N[0,1,2,3]
    // example for N = 12, partsize=4, size=3, rank=1, operates on N[(1*4)+0, (1*4)+1, (1*4)+2, (1*4)+3] --> N[4,5,6,7]
    // example for N = 12, partsize=4, size=3, rank=2, operates on N[(2*4)+0, (2*4)+1, (2*4)+2, (2*4)+3] --> N[8,9,10,11]

    // printf("rank %d working on array %lld exusive %lld\n", rank, start, stop);

    //finish send operations before entering loop
    if (start-1 >= 0) {
      // printf("send from %d to %d\n", rank, rank-1);
      MPI_Send(&(A[start]), 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
      }
    if (stop +1 < N){
      // printf("send from %d to %d\n", rank, rank+1);
      MPI_Send(&(A[stop-1]), 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
    }
    for (long long i = start; i < stop; i++) {
      //receive values before operating
      if (i == start && i-1 >= 0) {
        // printf("receive from %d to me=%d \n", rank-1, rank, i-1);
        MPI_Recv(&(A[i-1]), 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      if (i == stop -1 && i+1 < N){
        // printf("receive from %d to me=%d for node %d\n", rank+1, rank, i+1);
        MPI_Recv(&(A[i+1]), 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      // center stays constant (the heat is still on)
      if (i == source_x) {
        B[i] = A[i];
        continue;
      }

      // get temperature at current position
      value_t tc = A[i];

      // get temperatures of adjacent cells
      value_t tl = (i != 0) ? A[i - 1] : tc;
      value_t tr = (i != N - 1) ? A[i + 1] : tc;

      // compute new temperature at current position
      B[i] = tc + 0.2 * (tl + tr + (-2 * tc));

    }

    // swap matrices (just pointers, not content)
    Vector H = A;
    A = B;
    B = H;

    // show intermediate step
    if (!(t % 1000)) {
      printf("Step t=%d:\t", t);
      printTemperature(A, N);
      printf("\n");
    }
  }

  releaseVector(B);

  MPI_Gather(&A[start], partsize, MPI_DOUBLE, A, partsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // ---------- check ----------

  int success = 1;
  if (rank == 0){
    printf("Final:\t\t");
    printTemperature(A, N);
    printf("\n");
    for (long long i = 0; i < N; i++) {
      value_t temp = A[i];
      if (273 <= temp && temp <= 273 + 60)
        continue;
      success = 0;
      break;
    }
    printf("Verification: %s\n", (success) ? "OK" : "FAILED");
  }

  // ---------- cleanup ----------

  releaseVector(A);

  // done
  MPI_Finalize();
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Vector createVector(int N) {
  // create data and index vector
  return malloc(sizeof(value_t) * N);
}

void releaseVector(Vector m) { free(m); }

void printTemperature(Vector m, int N) {
  const char *colors = " .-:=+*^X#%@";
  const int numColors = 12;

  // boundaries for temperature (for simplicity hard-coded)
  const value_t max = 273 + 30;
  const value_t min = 273 + 0;

  // set the 'render' resolution
  int W = RESOLUTION;

  // step size in each dimension
  int sW = N / W;

  // room
  // left wall
  printf("X");
  // actual room
  for (int i = 0; i < W; i++) {
    // get max temperature in this tile
    value_t max_t = 0;
    for (int x = sW * i; x < sW * i + sW; x++) {
      max_t = (max_t < m[x]) ? m[x] : max_t;
    }
    value_t temp = max_t;

    // pick the 'color'
    int c = ((temp - min) / (max - min)) * numColors;
    c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);

    // print the average temperature
    printf("%c", colors[c]);
  }
  // right wall
  printf("X");
}
