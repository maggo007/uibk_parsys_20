CC=gcc
CC_FLAGS=-O2 -g -std=c99 -Wall -Wextra -pedantic

.PHONEY: all
all: heat_stencil_2D_seq_mpi heat_stencil_2D_omp

heat_stencil_2D_seq: heat_stencil_2D_seq.c
	@$(CC) $(CC_FLAGS) $< -o heat_stencil_2D_seq

# heat_stencil_2D_mpi: heat_stencil_2D_mpi.c
# 	@$(CC) $(CC_FLAGS) $< -o heat_stencil_2D_mpi

heat_stencil_2D_omp: heat_stencil_2D_omp.c
	@$(CC) $(CC_FLAGS) $< -o heat_stencil_2D_omp -fopenmp


.PHONEY: clean
clean:
	@rm heat_stencil_2D_seq
	@rm heat_stencil_2D_mpi
	@rm heat_stencil_2D_omp

.PHONEY: run
run: all
	@echo "Sequential:"
	@./heat_stencil_2D_seq
	@echo
