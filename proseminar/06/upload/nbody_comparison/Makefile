CC=gcc
CC_FLAGS=-O2 -std=gnu99 -Wall -Wextra -pedantic -lm -g -ffast-math

.PHONEY: all
all: nbody_seq nbody_omp nbody_omp_race nbody_bhut

nbody_seq: nbody_seq.c
	@$(CC) $(CC_FLAGS) $< -o nbody_seq

nbody_omp: nbody_omp.c
	@$(CC) $(CC_FLAGS) $< -o nbody_omp -fopenmp

nbody_omp_race: nbody_omp_race.c
	@$(CC) $(CC_FLAGS) $< -o nbody_omp_race -fopenmp

nbody_bhut: nbody_bhut.c
	@$(CC) $(CC_FLAGS) $< -o nbody_bhut

.PHONEY: clean
clean:
	@rm nbody_seq
	@rm nbody_omp
	@rm nbody_omp_race
	@rm nbody_bhut

.PHONEY: run
run: all
	@echo "seq:"
	@./nbody_seq
	@echo
