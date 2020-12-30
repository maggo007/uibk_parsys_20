CC=mpicc
CC_FLAGS=-O2 -g -std=gnu99 -Wall -Wextra -pedantic

.PHONEY: all
all: heat_stencil_2D_seq heat_stencil_2D_omp heat_stencil_2D_mpi heat_stencil_2D_mpi_hybrid heat_stencil_2D_mpi_shared

heat_stencil_2D_seq: heat_stencil_2D_seq.c
	@$(CC) $(CC_FLAGS) $< -o heat_stencil_2D_seq

 heat_stencil_2D_mpi: heat_stencil_2D_mpi.c
	@$(CC) $(CC_FLAGS) $< -o heat_stencil_2D_mpi

heat_stencil_2D_omp: heat_stencil_2D_omp.c
	@$(CC) $(CC_FLAGS) -fopenmp $< -o heat_stencil_2D_omp

heat_stencil_2D_mpi_hybrid: heat_stencil_2D_mpi_hybrid.c 
	@$(CC) $(CC_FLAGS) -fopenmp $< -o heat_stencil_2D_mpi_hybrid

heat_stencil_2D_mpi_shared: heat_stencil_2D_mpi_shared.c 
	@$(CC) $(CC_FLAGS) $< -o heat_stencil_2D_mpi_shared

.PHONEY: clean
clean:
	@rm heat_stencil_2D_seq
	@rm heat_stencil_2D_mpi
	@rm heat_stencil_2D_omp
	@rm heat_stencil_2D_mpi_hybrid
	@rm heat_stencil_2D_mpi_shared
