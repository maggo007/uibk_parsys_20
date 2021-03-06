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

  // blocking approach //

#pragma omp parallel for collapse(2)
  for (long long i = 0; i < N; i++) {
    for (long long j = 0; j < N; j++) {
      value_t sum = 0;
      for (long long k = 0; k < N; k++) {
        sum += A[i * N + k] * Btranspose[j * N + k];
      }
      C[i * N + j] = sum;
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
