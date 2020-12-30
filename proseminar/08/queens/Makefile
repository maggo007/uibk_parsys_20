CC=gcc
CC_FLAGS=-O2 -g -std=gnu99 -Wall -Wextra -pedantic

.PHONEY: all
all: queens_seq

queens_seq: queens_seq.c
	@$(CC) $(CC_FLAGS) -fopenmp $< -o queens_seq

.PHONEY: clean
clean:
	@rm queens_seq
