CC=gcc
CC_FLAGS=-O3 -g -pg -Wall -Wextra -pedantic -std=gnu99

.PHONEY: all
all: real real_omp

SRC=$(wildcard *.c)

real: $(SRC)
	@$(CC) $(CC_FLAGS) $(SRC) -o real -lm

real_omp: $(SRC)
	@$(CC) $(CC_FLAGS) -fopenmp $(SRC) -o real_omp -lm

.PHONEY: clean
clean:
	@rm real
	@rm real_omp

.PHONEY: run
run: all
	@echo "Sequential:"
	@perf stat -d ./real
	@echo "----------------------"
	@echo "Parallel:"
	@perf stat -d ./real_omp
	@echo
