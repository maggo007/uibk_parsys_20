//https://github.com/Spartee/Matrix-Multiplication-Chapel/blob/master/parmMultiply.chpl

//using CHPL_RT_NUM_THREADS_PER_LOCALE=  for number of threads. not working with 1

use Random;
use LinearAlgebra;

config const N = 2552;
config const seed = 314159265;
config const tasks = here.maxTaskPar;

//
// Output simulation setup.
//
writeln("Number of locales   = ", numLocales);
writeln("Number of points    = ", n);
writeln("Random number seed  = ", seed);
writeln("Number of tasks     = ", tasks, " (per locale)");

var A : [1..N, 1..N] real;  // first matrix
fillRandom(A, seed);

var B : [1..N, 1..N] real; // second matrix
forall (i,j) in B.domain do
  if (i==j) then B(i,j) = 1.0;

var C : [1..N, 1..N] real; // matrix for results.


coforall (i,j) in A.domain do {
  forall k in 1..N do {
    C[i,j] += A[i,k] * B[k,j];
  }
}
