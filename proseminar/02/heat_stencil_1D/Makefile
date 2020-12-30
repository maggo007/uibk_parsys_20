CC=mpicc
CC_FLAGS=-O2 -std=c99 -Wall -Wextra -pedantic

.PHONEY: all
all: heat_stencil_1D_seq heat_stencil_1D_mpi heat_stencil_1D_mpi_waves heat_stencil_1D_mpi_optimized

heat_stencil_1D_seq: heat_stencil_1D_seq.c
	@$(CC) $(CC_FLAGS) $< -o heat_stencil_1D_seq

heat_stencil_1D_mpi: heat_stencil_1D_mpi.c
	@$(CC) $(CC_FLAGS) $< -o heat_stencil_1D_mpi

heat_stencil_1D_mpi_waves: heat_stencil_1D_mpi_waves.c
	@$(CC) $(CC_FLAGS) $< -o heat_stencil_1D_mpi_waves

heat_stencil_1D_mpi_optimized: heat_stencil_1D_mpi_optimized.c
	@$(CC) $(CC_FLAGS) $< -o heat_stencil_1D_mpi_optimized

.PHONEY: clean
clean:
	@rm heat_stencil_1D_seq
	@rm heat_stencil_1D_mpi
	@rm heat_stencil_1D_mpi_waves
	@rm heat_stencil_1D_mpi_optimized

.PHONEY: run
run: all
	@echo "Sequential:"
	@./heat_stencil_1D_seq
	@echo
