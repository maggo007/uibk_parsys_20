//some inspritation from https://www.geeksforgeeks.org/n-queen-problem-backtracking-3/

#include <stdio.h>
#include <stdlib.h>

int check_paths(int pos_row, int pos_col, int problem_size, int board[problem_size][problem_size]){

  int row, col;
  //check left row
  for (int col = 0; col < pos_col; ++col) {
    if (board[pos_row][col]){
      //printf("left row\n");
      return 0;
    }
  }

  //check left upper diogonal
  for (row=pos_row-1, col=pos_col-1; row >= 0 && col >=0; --row, --col) {
    if (board[row][col]){
      //printf("left upper diogonal\n");
      return 0;
    }
  }

  //check left lower diogonal
  for (row=pos_row+1, col=pos_col-1; row < problem_size && col >=0; ++row, --col) {
    if (board[row][col]){
      //printf("left lower diogonal\n");
      return 0;
    }
  }
  return 1;
}

int req_solve(int problem_size, int board[problem_size][problem_size], int start_col, int lastrow, int *solvenumbers){
  if (start_col >= problem_size){
    /* for (int row=0; row < problem_size; ++row) { */
    /*   for (int col=0; col < problem_size; ++col) { */
    /*     printf("%d ", board[row][col]); */
    /*   } */
    /*   printf("\n"); */
    /* } */
    /* printf("Solution %d\n", ++(*solvenumbers)); */
    /* printf("---------\n"); */
#pragma omp atomic
    ++(*solvenumbers);
    return 0; //stop here normaly with return 1 for he first solution only
  }

  for (int row=0; row<problem_size; ++row) {
    if (row == lastrow || row == lastrow-1 || row == lastrow +1){
      continue;
    }
    if (check_paths(row, start_col, problem_size, board)){
      board[row][start_col] = 1;
      req_solve(problem_size, board, start_col+1, row, solvenumbers);
      board[row][start_col] = 0;
    }
  }
  return 0;
}

int main(int argc, char *argv[])
{
  // 'parsing' optional input parameter = problem size
  int N = 8;
  if (argc > 1) {
    N = atoi(argv[1]);
  }
  printf("Computing queens problem with N=%d x %d\n", N,N);

  // board is used with [row][columns]
  int *boards[N];
  //int board[N][N];
  int solvenumbers = 0;

for (int i=0; i < N; ++i) {
  boards[i] = malloc(sizeof(int)*N*N);
  for (int row=0; row < N; ++row) {
    for (int col=0; col < N; ++col) {
      *(boards[i]+row*N+col) = 0;
    }
    *(boards[i]+i*N) = 1;
  }
}

/* for (int row=0; row < N; ++row) { */
/*   for (int col=0; col < N; ++col) { */
/*     printf("%d ", *(boards[1]+row*N+col)); */
/*   } */
/*   printf("\n"); */
/* } */


#pragma omp parallel for
for (int i = 0; i < N; ++i) {
  req_solve(N, boards[i], 1, i, &solvenumbers);
}
  
  printf("found solutions=%d\n", solvenumbers);
  
  /* for (int row=0; row < N; ++row) { */
  /*   for (int col=0; col < N; ++col) { */
  /*     printf("%d ", board[row][col]); */
  /*   } */
  /*   printf("\n"); */
  /* } */
  for (int i = 0; i < N; ++i) {
    free(boards[i]);
  }

    
  return EXIT_SUCCESS;
}

