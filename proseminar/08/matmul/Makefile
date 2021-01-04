CC=gcc
CC_FLAGS=-O3 -g -std=gnu99 -Wall -Wextra -pedantic -ffast-math

.PHONEY: all
all: matmul_omp matmul_seq matmul_omp_op matmul_omp_op_block

matmul_omp: matmul_omp.c
	@$(CC) $(CC_FLAGS) -fopenmp $< -o matmul_omp

matmul_seq: matmul_omp.c
	@$(CC) $(CC_FLAGS) $< -o matmul_seq

matmul_omp_op: matmul_omp_op.c
	@$(CC) $(CC_FLAGS) -fopenmp $< -o matmul_omp_op

matmul_omp_op_block: matmul_omp_op_block.c
	@$(CC) $(CC_FLAGS) -fopenmp $< -o matmul_omp_op_block

.PHONEY: clean
clean:
	@rm matmul_omp
	@rm matmul_seq
	@rm matmul_omp_op
	@rm matmul_omp_op_block
