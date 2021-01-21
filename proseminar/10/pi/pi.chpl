//https://github.com/chapel-lang/chapel/blob/master/test/exercises/MonteCarloPi/deitz/MonteCarloPi

use Random;

config const problemSize = 1000000000;
//config const problemSize = 10000000;
config const seed = 314159265;

writeln("Number of points    = ", problemSize);
writeln("Random number seed  = ", seed);

var rs = new owned NPBRandomStream(real, seed, parSafe=false);
var count = 0;
for i in 1..problemSize do
  count += rs.getNext()**2 + rs.getNext()**2 <= 1.0;

writef("Approximation of pi = %{#.#######}\n", count * 4.0 / problemSize);
