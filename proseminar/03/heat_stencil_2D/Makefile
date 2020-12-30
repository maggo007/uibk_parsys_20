CC=mpicc
CC_FLAGS=-O2 -g -std=c99 -Wall -Wextra -pedantic

.PHONEY: all
all: heat_stencil_2D_seq heat_stencil_2D_mpi heat_stencil_2D_mpi_optimized2

heat_stencil_2D_seq: heat_stencil_2D_seq.c
	@$(CC) $(CC_FLAGS) $< -o heat_stencil_2D_seq

heat_stencil_2D_mpi: heat_stencil_2D_mpi.c
	@$(CC) $(CC_FLAGS) $< -o heat_stencil_2D_mpi

# heat_stencil_2D_mpi_optimized: heat_stencil_2D_mpi_optimized.c
# 	@$(CC) $(CC_FLAGS) $< -o heat_stencil_2D_mpi_optimized

heat_stencil_2D_mpi_optimized2: heat_stencil_2D_mpi_optimized2.c
	@$(CC) $(CC_FLAGS) $< -o heat_stencil_2D_mpi_optimized2

.PHONEY: clean
clean:
	@rm heat_stencil_2D_seq
	@rm heat_stencil_2D_mpi
	@rm heat_stencil_2D_mpi_optimized
	@rm heat_stencil_2D_mpi_optimized2

.PHONEY: run
run: all
	@echo "Sequential:"
	@./heat_stencil_2D_seq
	@echo
