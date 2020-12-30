CC=gcc
CC_FLAGS=-O2 -std=gnu99 -Wall -Wextra -pedantic

.PHONEY: all
all: pi_omp pi_seq pi_omp_reduction pi_omp_atomic pi_omp_reduction2 pi_omp_critical pi_omp_reduction2_float pi_omp_reduction_float pi_omp_atomic_float pi_omp_critical_float

pi_seq: pi_seq.c
	@$(CC) $(CC_FLAGS) $< -o pi_seq

pi_omp: pi_omp.c
	@$(CC) $(CC_FLAGS) $< -o pi_omp -fopenmp

pi_omp_reduction: pi_omp_reduction.c
	@$(CC) $(CC_FLAGS) $< -o pi_omp_reduction -fopenmp

pi_omp_reduction_float: pi_omp_reduction_float.c
	@$(CC) $(CC_FLAGS) $< -o pi_omp_reduction_float -fopenmp

pi_omp_reduction2: pi_omp_reduction2.c
	@$(CC) $(CC_FLAGS) $< -o pi_omp_reduction2 -fopenmp

pi_omp_reduction2_float: pi_omp_reduction2_float.c
	@$(CC) $(CC_FLAGS) $< -o pi_omp_reduction2_float -fopenmp

pi_omp_atomic: pi_omp_atomic.c
	@$(CC) $(CC_FLAGS) $< -o pi_omp_atomic -fopenmp

pi_omp_atomic_float: pi_omp_atomic_float.c
	@$(CC) $(CC_FLAGS) $< -o pi_omp_atomic_float -fopenmp

pi_omp_critical: pi_omp_critical.c
	@$(CC) $(CC_FLAGS) $< -o pi_omp_critical -fopenmp

pi_omp_critical_float: pi_omp_critical_float.c
	@$(CC) $(CC_FLAGS) $< -o pi_omp_critical_float -fopenmp

.PHONEY: clean
clean:
	@rm pi_omp
	@rm pi_seq
	@rm pi_omp_atomic
	@rm pi_omp_reduction
	@rm pi_omp_reduction2
	@rm pi_omp_critical
	@rm pi_omp_reduction_float
	@rm pi_omp_reduction2_float
	@rm pi_omp_atomic_float
	@rm pi_omp_critical_float

.PHONEY: run
run: all
	@echo "openmp:"
	@./pi_omp
	@echo
