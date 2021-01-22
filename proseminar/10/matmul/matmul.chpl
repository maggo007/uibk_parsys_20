//https://github.com/Spartee/Matrix-Multiplication-Chapel/blob/master/parmMultiply.chpl


use Random;
use LinearAlgebra;

config const N = 2552;
config const seed = 314159265;

var A : [1..N, 1..N] real;  // first matrix
fillRandom(A, seed);

var B : [1..N, 1..N] real; // second matrix
forall (i,j) in B.domain do
  if (i==j) then B(i,j) = 1.0;

var C : [1..N, 1..N] real; // matrix for results.

//transposed access
forall (i,j) in A.domain do {
  forall k in 1..N do {
    C[i,j] += A[i,k] * B[j,k];
  }
}
