//some parts are copied from our old solutions in parallel openCL course https://git.uibk.ac.at/csat2062/parallel_local

#include <stdio.h>
#include <stdlib.h>

typedef double value_t;

#define RESOLUTION 120

// -- vector utilities --

typedef value_t *Vector;

Vector createVector(int N, int M);

void releaseVector(Vector m);

void printTemperature(Vector m, int N, int M);

// -- simulation code ---

int main(int argc, char **argv) {
  // 'parsing' optional input parameter = problem size
  int N = 1000;
  int M = 1000;
  if (argc == 3) {
    N = atoi(argv[1]);
    M = atoi(argv[2]);
  }
  int T = N * 500;
  printf("Computing heat-distribution for room size N=%d x %d for T=%d timesteps\n", N,M, T);

  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Vector A = createVector(N,M);

  // set up initial conditions in A
  for (int y = 0; y < M; ++y) {
    for (int x = 0; x < N; ++x) {
      A[x + y*N] = 273; // temperature is 0° C everywhere (273 K)
    }
  }

  // and there is a heat source in one corner
  int source_x = N / 4;
  int source_y = M / 4;
  A[source_y * N + source_x] = 273 + 60;

  printf("Initial:\n");
  printTemperature(A, N, M);
  printf("\n");

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(N,M);

  // for each time step ..
  for (int t = 0; t < T; t++) {
    // .. we propagate the temperature
    for (int y = 0; y < M; y++) {
      for (int x = 0; x < N; x++) {
        // center stays constant (the heat is still on)
        if (x == source_x && y == source_y) {
          B[y*N+x] = A[y*N+x];
          continue;
        }
        
        // get current temperature at (x,y)
        value_t tc = A[y*N+x];

        // get temperatures left/right and up/down
        value_t tl = ( x !=  0  ) ? A[y*N+(x-1)] : tc;
        value_t tr = ( x != N-1 ) ? A[y*N+(x+1)] : tc;
        value_t tu = ( y !=  0  ) ? A[(y-1)*N+x] : tc;
        value_t td = ( y != M-1 ) ? A[(y+1)*N+x] : tc;

        // update temperature at current point
        B[y*N+x] = tc + 1.0/5 * (tl + tr + tu + td + (-4*tc));

      }
    }

    // swap matrices (just pointers, not content)
    Vector H = A;
    A = B;
    B = H;

    // show intermediate step
    if (!(t % 1000)) {
      printf("Step t=%d:\n", t);
      printTemperature(A, N, M);
      printf("\n");
    }
  }

  releaseVector(B);

  // ---------- check ----------

  printf("Final:\n");
  printTemperature(A, N, M);
  printf("\n");

  int success = 1;
  for(int y = 0; y<M; y++) {
    for(int x = 0; x<N; x++) {
      value_t temp = A[y*N+x];
      if (273 <= temp && temp <= 273+60) continue;
      success = 0;
      break;
    }
  }


  printf("Verification: %s\n", (success) ? "OK" : "FAILED");

  // ---------- cleanup ----------

  releaseVector(A);

  // done
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Vector createVector(int N, int M) {
  // create data and index vector
  return malloc(sizeof(value_t) * N * M);
}

void releaseVector(Vector m) { free(m); }


//taken from old parallel course https://git.uibk.ac.at/csat2062/parallel_local
void printTemperature(Vector m, int N, int M) {
    const char* colors = " .-:=+*#%@";
    const int numColors = 10;

    // boundaries for temperature (for simplicity hard-coded)
    const value_t max = 273 + 30;
    const value_t min = 273 + 0;

    // set the 'render' resolution
    int H = 30;
    int W = 60;

    // step size in each dimension
    int sH = M/H;
    int sW = N/W;


    // upper wall
    for(int i=0; i<W+2; i++) {
        printf("X");
    }
    printf("\n");

    // room
    for(int i=0; i<H; i++) {
        // left wall
        printf("X");
        // actual room
        for(int j=0; j<W; j++) {

            // get max temperature in this tile
            value_t max_t = 0;
            for(int x=sH*i; x<sH*i+sH; x++) {
                for(int y=sW*j; y<sW*j+sW; y++) {
                    max_t = (max_t < m[x*N+y]) ? m[x*N+y] : max_t;
                }
            }
            value_t temp = max_t;

            // pick the 'color'
            int c = ((temp - min) / (max - min)) * numColors;
            c = (c >= numColors) ? numColors-1 : ((c < 0) ? 0 : c);

            // print the average temperature
            printf("%c",colors[c]);
        }
        // right wall
        printf("X\n");
    }

    // lower wall
    for(int i=0; i<W+2; i++) {
        printf("X");
    }
    printf("\n");

}
