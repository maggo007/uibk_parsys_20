// some parts are copied from our old solutions in parallel openCL course
// https://git.uibk.ac.at/csat2062/parallel_local

#include <stdio.h>
#include <stdlib.h>

typedef double value_t;

#define RESOLUTION 120
#define NDEBUG

// -- vector utilities --

typedef value_t *Vector;

Vector createVector(int N, int M, int O);

void releaseVector(Vector m);

void printTemperature(Vector m, int N, int M, int O);

// -- simulation code ---

int main(int argc, char **argv) {
  // 'parsing' optional input parameter = problem size
  int N = 100;
  int M = 100;
  int Oz = 100;
  if (argc >1) {
    N = atoi(argv[1]);
    M = atoi(argv[1]);
    Oz = atoi(argv[1]);
  }
  int T = 50000;
  printf("Computing heat-distribution for room size N=%d x %d x %d for T=%d timesteps\n", N, M, Oz, T);

  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Vector A = createVector(N, M, Oz);

  // set up initial conditions in A
  for (int z = 0; z < Oz; ++z) {
    for (int y = 0; y < M; ++y) {
      for (int x = 0; x < N; ++x) {
        A[x + y * N + z * N * M] =  273; // temperature is 0° C everywhere (273 K)
      }
    }
  }

  // and there is a heat source in one corner
  int source_x = N / 4;
  int source_y = M / 4;
  int source_z = Oz / 4;
  A[source_z*N*M + source_y*N + source_x] = 273 + 60;

  printf("Initial:\n");
  printTemperature(A, N, M, Oz);
  printf("\n");

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(N, M, Oz);

  // for each time step ..
  for (int t = 0; t < T; t++) {
    // .. we propagate the temperature
    for (int z = 0; z < Oz; z++) {
      for (int y = 0; y < M; y++) {
        for (int x = 0; x < N; x++) {
          // center stays constant (the heat is still on)
          if (x == source_x && y == source_y && z == source_z) {
            B[z*N*M + y*N + x] = A[z*N*M + y*N + x];
#ifndef NDEBUG
            printf("written origin %f on x=%d y=%d z=%d on step %d\n",A[x + y * N + z * N * M], x, y, z, t);
#endif
          
            continue;
          }

          // get current temperature at (x,y,z)
          value_t tc = A[z*N*M + y*N + x];

          // get temperatures left/right and up/down  and below/above
          value_t tl = (x != 0    ) ? A[z*N*M + y*N + (x-1)] : tc;
          value_t tr = (x != N - 1) ? A[z*N*M + y*N + (x+1)] : tc;
          value_t tu = (y != 0    ) ? A[z*N*M + (y-1)*N + x] : tc;
          value_t td = (y != M - 1) ? A[z*N*M + (y+1)*N + x] : tc;
          value_t tb = (z != 0    ) ? A[(z-1)*N*M + y*N + x] : tc;
          value_t ta = (z != Oz - 1) ? A[(z+1)*N*M + y*N + x] : tc;

          // update temperature at current point
          double temp = tc + 1.0/7 * (tl + tr + tu + td + tb + ta + (-6 * tc));
#ifndef NDEBUG
          if (tc<0 || temp < 0){
            printf("negative value\n");
          }
          if (temp > 273+60){
            printf("to high value\n");
          }
          if (temp < 273){
            printf("below start value\n");
          }
          printf("tc=%f\n", tc);
          printf("tl=%f\n", tl);
          printf("tr=%f\n", tr);
          printf("tu=%f\n", tu);
          printf("td=%f\n", td);
          printf("tb=%f\n", tb);
          printf("ta=%f\n", ta);
          printf("written %f on x=%d y=%d z=%d on step %d\n", temp, x, y, z, t);
          printf("---\n");
#endif
          B[z*N*M + y*N + x] = temp;
        }
      }
    }
    // swap matrices (just pointers, not content)
    Vector H = A;
    A = B;
    B = H;

    // show intermediate step
    if (!(t % 1000) && 0) {
      printf("Step t=%d:\n", t);
      printTemperature(A, N, M, Oz);
      printf("\n");
    }
  }

  releaseVector(B);

  // ---------- check ----------

  printf("Final:\n");
  printTemperature(A, N, M, Oz);
  printf("\n");

  int success = 1;
  for (int z = 0; z < Oz; z++) {
    for (int y = 0; y < M; y++) {
      for (int x = 0; x < N; x++) {
        value_t temp = A[z*N*M + y*N + x];
        if (273 <= temp && temp <= 273 + 60)
          continue;
        success = 0;
        printf("failure on cell x=%d y=%d z=%d with value %f\n", x, y, z, temp);
        break;
      }
    }
  }

  printf("Verification: %s\n", (success) ? "OK" : "FAILED");

  // ---------- cleanup ----------

  releaseVector(A);

  // done
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Vector createVector(int N, int M, int O) {
  // create data and index vector
  return malloc(sizeof(value_t) * N * M * O);
}

void releaseVector(Vector m) { free(m); }

// taken from old parallel course https://git.uibk.ac.at/csat2062/parallel_local
void printTemperature(Vector m, int N, int M, int Oz) {
  const char *colors = " .-:=+*#%@";
  const int numColors = 10;

  // boundaries for temperature (for simplicity hard-coded)
  const value_t max = 273 + 30;
  const value_t min = 273 + 0;

  // set the 'render' resolution
  int H = 30;
  int W = 50;
  int D = 50;

  // step size in each dimension
  int sH = N / H;
  int sW = M / W;
  int sD = Oz / D;

  // upper wall
  for (int i = 0; i < W + 2; i++) {
    printf("X");
  }
  printf("\n");

  // room
  for (int i = 0; i < H; i++) {
    // left wall
    printf("X");
    // actual room
    for (int j = 0; j < W; j++) {

      // get max temperature in this tile
      value_t max_t = 0;
      for (int x = sH * i; x < sH * i + sH; x++) {
        for (int y = sW * j; y < sW * j + sW; y++) {
          //add values over z coordinate
          for (int z = 0; z < Oz; z++) {
            max_t = (max_t < m[x * N + y + N*M*z]) ? m[x * N + y+ N*M*z] : max_t;
          }
        }
      }
      // avarage over z hight
      value_t temp = max_t;

      // pick the 'color'
      int c = ((temp - min) / (max - min)) * numColors;
      c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);

      // print the average temperature
      printf("%c", colors[c]);
    }
    // right wall
    printf("X\n");
  }

  // lower wall
  for (int i = 0; i < W + 2; i++) {
    printf("X");
  }
  printf("\n");
}
