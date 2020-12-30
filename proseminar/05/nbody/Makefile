CC=gcc
CC_FLAGS=-O2 -std=gnu99 -Wall -Wextra -pedantic -lm -ffast-math

.PHONEY: all
all: nbody_seq nbody_seq2 nbody_omp

nbody_seq: nbody_seq.c
	@$(CC) $(CC_FLAGS) $< -o nbody_seq

nbody_seq2: nbody_seq2.c
	@$(CC) $(CC_FLAGS) $< -o nbody_seq2

nbody_omp: nbody_omp.c
	@$(CC) $(CC_FLAGS) $< -o nbody_omp -fopenmp

.PHONEY: clean
clean:
	@rm nbody_seq
	@rm nbody_omp

.PHONEY: run
run: all
	@echo "seq:"
	@./nbody_seq
	@echo
