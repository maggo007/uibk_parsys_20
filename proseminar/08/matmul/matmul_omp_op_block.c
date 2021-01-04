// copy from openCL course  https://git.uibk.ac.at/csat2062/parallel_local/

#include <stdio.h>
#include <stdlib.h>

typedef double value_t;

// -- matrix utilities --

typedef value_t *Matrix;

Matrix createMatrix(int N, int M);

void releaseMatrix(Matrix m);

// ----------------------

int main(int argc, char **argv) {

  // 'parsing' optional input parameter = problem size
  int N = 1000;
  if (argc > 1) {
    N = atoi(argv[1]);
  }
  printf("Computing matrix-matrix product with N=%d\n", N);

  // ---------- setup ----------

  // create two input matrixes (on heap!)
  Matrix A = createMatrix(N, N);
  Matrix B = createMatrix(N, N);
  Matrix Btranspose = createMatrix(N, N);

  // fill matrixes
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      A[i * N + j] = i * j;            // some matrix - note: flattend indexing!
      B[i * N + j] = (i == j) ? 1 : 0; // identity
    }
  }

  // transpose B matrix for better data coherence //
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      Btranspose[i * N + j] = B[j * N + i];
    }
  }
  

  // ---------- compute ----------

  Matrix C = createMatrix(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      C[i * N + j] = 0.0;
    }
  }

  // blocking approach //

  int blocksize = 60;
  value_t sum = 0;
#pragma omp parallel for collapse(2)
  for (long long ii = 0; ii < N; ii+=blocksize) {
    for (long long jj = 0; jj < N; jj+=blocksize) {
      for (long long kk = 0; kk < N; kk+=blocksize) {
        for (long long i = ii; i < ii+blocksize; i++) {
          for (long long j = jj; j < jj+blocksize; j++) {
            sum = C[i*N +j];
            for (long long k = kk; k < kk+blocksize; k++) {
              sum += A[i * N + k] * Btranspose[j * N + k];
            }
            C[i * N + j] = sum;
          }
        }
      }
    }
  }

  // ---------- check ----------

  int success = 1;
  for (long long i = 0; i < N; i++) {
    for (long long j = 0; j < N; j++) {
      if (C[i * N + j] == i * j)
        continue;
      success = 0;
      break;
    }
  }

  printf("Verification: %s\n", (success) ? "OK" : "FAILED");
  printf("blocksize=%d \n", blocksize);

  // ---------- cleanup ----------

  releaseMatrix(A);
  releaseMatrix(B);
  releaseMatrix(C);

  // done
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Matrix createMatrix(int N, int M) {
  // create data and index vector
  return malloc(sizeof(value_t) * N * M);
}

void releaseMatrix(Matrix m) { free(m); }
