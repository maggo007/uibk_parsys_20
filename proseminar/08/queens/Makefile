CC=gcc
CC_FLAGS=-O3 -g -std=gnu99 -Wall -Wextra -pedantic

.PHONEY: all
all: queens_seq queens_omp queens_omp_op

queens_seq: queens_seq.c
	@$(CC) $(CC_FLAGS) -fopenmp $< -o queens_seq

queens_omp: queens_omp.c
	@$(CC) $(CC_FLAGS) -fopenmp $< -o queens_omp

queens_omp_op: queens_omp_op.c
	@$(CC) $(CC_FLAGS) -fopenmp $< -o queens_omp_op

.PHONEY: clean
clean:
	@rm queens_seq
	@rm queens_omp
	@rm queens_omp_op
