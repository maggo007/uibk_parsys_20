CC=gcc
CC_FLAGS=-O3 -g -std=gnu99 -Wall -Wextra -pedantic

.PHONEY: all
all: matmul_omp matmul_seq

matmul_omp: matmul_omp.c
	@$(CC) $(CC_FLAGS) -fopenmp $< -o matmul_omp

matmul_seq: matmul_omp.c
	@$(CC) $(CC_FLAGS) $< -o matmul_seq

.PHONEY: clean
clean:
	@rm matmul_omp
	@rm matmul_seq
