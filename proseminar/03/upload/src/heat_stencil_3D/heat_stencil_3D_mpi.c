//some parts are copied from our old solutions in parallel openCL course https://git.uibk.ac.at/csat2062/parallel_local

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

typedef double value_t;

#define RESOLUTION 120
#define DIMENSIONS 3

enum mode {SEND, RECEIVE};

// -- vector utilities --

typedef value_t *Vector;

Vector createVector(int N, int M, int Oz);

void releaseVector(Vector m);

void printTemperature(Vector m, int N, int M, int Oz);

void neighbourMPI(int dimension, int* arraysize, int* subsize, int* substart, int* neighbour, MPI_Comm cartesian, double* sendbuffer, enum mode mode);

// -- simulation code ---
int periodic[DIMENSIONS], coord[DIMENSIONS];

int main(int argc, char **argv) {
  int rank;
  int size;
  int partsize_row;
  int partsize_col;
  int partsize_slide;
  int dim[DIMENSIONS] = {0,0,0};
  
  MPI_Comm cart;
  //MPI_Datatype subarray;
  MPI_Datatype gather, gather2;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  printf("rank=%d with size=%d\n", rank, size);
  // 'parsing' optional input parameter = problem size
  int N = 300;
  if (argc > 1) {
    N = atoi(argv[1]);
  }
  int T = 50000;
  
  /* if (N % size != 0) { */
  /*   printf("size not ok for problem size. not dividable without remainder\n"); */
  /*   exit(EXIT_FAILURE); */
  /* } */

  int arraysize[DIMENSIONS] = {N,N,N};

  printf("Computing heat-distribution for room size N=%d x %d x %d for T=%d timesteps\n", N,N,N, T);

  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Vector A = createVector(N,N,N);

  // set up initial conditions in A
  for (int slide = 0; slide < N; ++slide) {
    for (int row = 0; row < N; ++row) {
      for (int col = 0; col < N; ++col) {
        A[col + row*N + slide*N*N] = 273; // temperature is 0° C everywhere (273 K)
      }
    }
  }

  MPI_Dims_create(size, DIMENSIONS, dim);
  printf("dims=%d %d %d\n", dim[0], dim[1], dim[2]);
  if (N % dim[0] != 0 || N % dim[1] != 0 || N % dim[2] != 0){
    printf("size not ok for problem size. not dividable without remainder\n");
    exit(EXIT_FAILURE);    
  }
  partsize_row = N / dim[0];
  partsize_col = N / dim[1];
  partsize_slide = N / dim[2];
  //dim[0] = size/2;
  //dim[1] = size/2;
  periodic[0] = 0;
  periodic[1] = 0;
  periodic[2] = 0;
  int reorder = 0;
  MPI_Cart_create(MPI_COMM_WORLD, DIMENSIONS, dim, periodic, reorder, &cart);

  MPI_Cart_coords(cart, rank, DIMENSIONS, coord);

  // and there is a heat source in one corner
  int source_row = N / 4;
  int source_col = N / 4;
  int source_slide = N / 4;
  A[source_row * N + source_col + N*N* source_slide] = 273 + 60;

  printf("Initial:\n");
  printTemperature(A, N, N,N);
  printf("\n");

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(N,N,N);

  for (int slide = 0; slide < N; ++slide) {
    for (int row = 0; row < N; ++row) {
      for (int col = 0; col < N; ++col) {
        B[col + row*N + slide*N*N] = 0; // temperature is 0° C everywhere (273 K)
      }
    }
  }

  long long start_row = coord[0] * partsize_row;
  long long stop_row = coord[0] * partsize_row + partsize_row;
  long long start_col = coord[1] * partsize_col;
  long long stop_col = coord[1] * partsize_col + partsize_col;
  long long start_slide = coord[2] * partsize_slide;
  long long stop_slide = coord[2] * partsize_slide + partsize_slide;

  printf("rank %d is at coord(%d,%d, %d) with array start_row=%lld stop_row=%lld start_col=%lld stop_col=%lld start_slide=%lld stop_slide=%lld \n", rank, coord[0], coord[1], coord[2], start_row, stop_row, start_col, stop_col, start_slide, stop_slide);
  int neighbour[DIMENSIONS];
  int subsize[DIMENSIONS];
  int substart[DIMENSIONS];
  // for each time step ..
  for (int t = 0; t < T; t++) {
    // .. we propagate the temperature
        
    if((coord[0]+coord[1]+coord[2]) % 2 == 0){
      // send to left neighbour
      if (coord[1] != 0){
        subsize[0] = partsize_row;
        subsize[1] = 1;
        subsize[2] = partsize_slide;
        substart[0] = start_row;
        substart[1] = start_col;
        substart[2] = start_slide;
        neighbour[0] = coord[0];
        neighbour[1] = coord[1]-1;
        neighbour[2] = coord[2];
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, SEND);
      }

      // receive from right neighbour
      if (coord[1] < dim[1]-1){
        subsize[0] = partsize_row;
        subsize[1] = 1;
        subsize[2] = partsize_slide;
        substart[0] = start_row;
        substart[1] = stop_col;
        substart[2] = start_slide;
        neighbour[0] = coord[0];
        neighbour[1] = coord[1]+1;
        neighbour[2] = coord[2];
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, RECEIVE);
      }
      
      // send to right  neighbour
      if (coord[1] < dim[1]-1){
        subsize[0] = partsize_row;
        subsize[1] = 1;
        subsize[2] = partsize_slide;
        substart[0] = start_row;
        substart[1] = stop_col-1;
        substart[2] = start_slide;
        neighbour[0] = coord[0];
        neighbour[1] = coord[1]+1;
        neighbour[2] = coord[2];
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, SEND);
      }

      // receive from left neighbour
      if (coord[1] != 0){
        subsize[0] = partsize_row;
        subsize[1] = 1;
        subsize[2] = partsize_slide;
        substart[0] = start_row;
        substart[1] = start_col-1;
        substart[2] = start_slide;
        neighbour[0] = coord[0];
        neighbour[1] = coord[1]-1;
        neighbour[2] = coord[2];
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, RECEIVE);
      }
      
      // send to bottom neighbour
      if (coord[0] < dim[0]-1){
        subsize[0] = 1;
        subsize[1] = partsize_col;
        subsize[2] = partsize_slide;
        substart[0] = stop_row-1;
        substart[1] = start_col;
        substart[2] = start_slide;
        neighbour[0] = coord[0]+1;
        neighbour[1] = coord[1];
        neighbour[2] = coord[2];
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, SEND);
      }
      // receive from top neighbour
      if (coord[0] != 0){
        subsize[0] = 1;
        subsize[1] = partsize_col;
        subsize[2] = partsize_slide;
        substart[0] = start_row-1;
        substart[1] = start_col;
        substart[2] = start_slide;
        neighbour[0] = coord[0]-1;
        neighbour[1] = coord[1];
        neighbour[2] = coord[2];
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, RECEIVE);
      }
      // send to top neighbour
      if (coord[0] != 0){
        subsize[0] = 1;
        subsize[1] = partsize_col;
        subsize[2] = partsize_slide;
        substart[0] = start_row;
        substart[1] = start_col;
        substart[2] = start_slide;
        neighbour[0] = coord[0]-1;
        neighbour[1] = coord[1];
        neighbour[2] = coord[2];
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, SEND);
      }
      // receive from bottom neighbour
      if (coord[0] < dim[0]-1){
        subsize[0] = 1;
        subsize[1] = partsize_col;
        subsize[2] = partsize_slide;
        substart[0] = stop_row;
        substart[1] = start_col;
        substart[2] = start_slide;
        neighbour[0] = coord[0]+1;
        neighbour[1] = coord[1];
        neighbour[2] = coord[2];
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, RECEIVE);
      }
      
      // send to below neighbour
      if (coord[2] != 0){
        subsize[0] = partsize_row;
        subsize[1] = partsize_col;
        subsize[2] = 1;
        substart[0] = start_row;
        substart[1] = start_col;
        substart[2] = start_slide;
        neighbour[0] = coord[0];
        neighbour[1] = coord[1];
        neighbour[2] = coord[2]-1;
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, SEND);
      }
      // receive from above neighbour
      if (coord[2] < dim[2]-1){
        subsize[0] = partsize_row;
        subsize[1] = partsize_col;
        subsize[2] = 1;
        substart[0] = start_row;
        substart[1] = start_col;
        substart[2] = stop_slide;
        neighbour[0] = coord[0];
        neighbour[1] = coord[1];
        neighbour[2] = coord[2]+1;
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, RECEIVE);
      }
      
      // send to above neighbour
      if (coord[2] < dim[2]-1){
        subsize[0] = partsize_row;
        subsize[1] = partsize_col;
        subsize[2] = 1;
        substart[0] = start_row;
        substart[1] = start_col;
        substart[2] = stop_slide-1;
        neighbour[0] = coord[0];
        neighbour[1] = coord[1];
        neighbour[2] = coord[2]+1;
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, SEND);
      }

      // receive from below neighbour
      if (coord[2] != 0 ){
        subsize[0] = partsize_row;
        subsize[1] = partsize_col;
        subsize[2] = 1;
        substart[0] = start_row;
        substart[1] = start_col;
        substart[2] = start_slide-1;
        neighbour[0] = coord[0];
        neighbour[1] = coord[1];
        neighbour[2] = coord[2]-1;
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, RECEIVE);
      }

    }else{

      // receive from right neighbour
      if (coord[1] < dim[1]-1){
        subsize[0] = partsize_row;
        subsize[1] = 1;
        subsize[2] = partsize_slide;
        substart[0] = start_row;
        substart[1] = stop_col;
        substart[2] = start_slide;
        neighbour[0] = coord[0];
        neighbour[1] = coord[1]+1;
        neighbour[2] = coord[2];
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, RECEIVE);
      }

      // send to left neighbour
      if (coord[1] != 0){
        subsize[0] = partsize_row;
        subsize[1] = 1;
        subsize[2] = partsize_slide;
        substart[0] = start_row;
        substart[1] = start_col;
        substart[2] = start_slide;
        neighbour[0] = coord[0];
        neighbour[1] = coord[1]-1;
        neighbour[2] = coord[2];
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, SEND);
      }
      
      
      // receive from left neighbour
      if (coord[1] != 0){
        subsize[0] = partsize_row;
        subsize[1] = 1;
        subsize[2] = partsize_slide;
        substart[0] = start_row;
        substart[1] = start_col-1;
        substart[2] = start_slide;
        neighbour[0] = coord[0];
        neighbour[1] = coord[1]-1;
        neighbour[2] = coord[2];
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, RECEIVE);
      }

      // send to right  neighbour
      if (coord[1] < dim[1]-1){
        subsize[0] = partsize_row;
        subsize[1] = 1;
        subsize[2] = partsize_slide;
        substart[0] = start_row;
        substart[1] = stop_col-1;
        substart[2] = start_slide;
        neighbour[0] = coord[0];
        neighbour[1] = coord[1]+1;
        neighbour[2] = coord[2];
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, SEND);
      }
      
      
      
      // receive from top neighbour
      if (coord[0] != 0){
        subsize[0] = 1;
        subsize[1] = partsize_col;
        subsize[2] = partsize_slide;
        substart[0] = start_row-1;
        substart[1] = start_col;
        substart[2] = start_slide;
        neighbour[0] = coord[0]-1;
        neighbour[1] = coord[1];
        neighbour[2] = coord[2];
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, RECEIVE);
      }
      // send to bottom neighbour
      if (coord[0] < dim[0]-1){
        subsize[0] = 1;
        subsize[1] = partsize_col;
        subsize[2] = partsize_slide;
        substart[0] = stop_row-1;
        substart[1] = start_col;
        substart[2] = start_slide;
        neighbour[0] = coord[0]+1;
        neighbour[1] = coord[1];
        neighbour[2] = coord[2];
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, SEND);
      }
      
      
      
      // receive from bottom neighbour
      if (coord[0] < dim[0]-1){
        subsize[0] = 1;
        subsize[1] = partsize_col;
        subsize[2] = partsize_slide;
        substart[0] = stop_row;
        substart[1] = start_col;
        substart[2] = start_slide;
        neighbour[0] = coord[0]+1;
        neighbour[1] = coord[1];
        neighbour[2] = coord[2];
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, RECEIVE);
      }
      // send to top neighbour
      if (coord[0] != 0){
        subsize[0] = 1;
        subsize[1] = partsize_col;
        subsize[2] = partsize_slide;
        substart[0] = start_row;
        substart[1] = start_col;
        substart[2] = start_slide;
        neighbour[0] = coord[0]-1;
        neighbour[1] = coord[1];
        neighbour[2] = coord[2];
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, SEND);
      }
      
      
      
      // receive from above neighbour
      if (coord[2] < dim[2]-1){
        subsize[0] = partsize_row;
        subsize[1] = partsize_col;
        subsize[2] = 1;
        substart[0] = start_row;
        substart[1] = start_col;
        substart[2] = stop_slide;
        neighbour[0] = coord[0];
        neighbour[1] = coord[1];
        neighbour[2] = coord[2]+1;
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, RECEIVE);
      }
      // send to below neighbour
      if (coord[2] != 0){
        subsize[0] = partsize_row;
        subsize[1] = partsize_col;
        subsize[2] = 1;
        substart[0] = start_row;
        substart[1] = start_col;
        substart[2] = start_slide;
        neighbour[0] = coord[0];
        neighbour[1] = coord[1];
        neighbour[2] = coord[2]-1;
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, SEND);
      }
      
      

      // receive from below neighbour
      if (coord[2] != 0 ){
        subsize[0] = partsize_row;
        subsize[1] = partsize_col;
        subsize[2] = 1;
        substart[0] = start_row;
        substart[1] = start_col;
        substart[2] = start_slide-1;
        neighbour[0] = coord[0];
        neighbour[1] = coord[1];
        neighbour[2] = coord[2]-1;
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, RECEIVE);
      }
      // send to above neighbour
      if (coord[2] < dim[2]-1){
        subsize[0] = partsize_row;
        subsize[1] = partsize_col;
        subsize[2] = 1;
        substart[0] = start_row;
        substart[1] = start_col;
        substart[2] = stop_slide-1;
        neighbour[0] = coord[0];
        neighbour[1] = coord[1];
        neighbour[2] = coord[2]+1;
        neighbourMPI(DIMENSIONS, arraysize, subsize, substart, neighbour, cart, A, SEND);
      }
      
    }
    
    //MPI_Barrier(MPI_COMM_WORLD);
    for (int slide = start_slide; slide < stop_slide; slide++){
      for (int row = start_row; row < stop_row; row++) {
        for (int col = start_col; col < stop_col; col++) {
          // center stays constant (the heat is still on)
          if (col == source_col && row == source_row && slide == source_slide) {
            B[slide*N*N+row*N+col] = A[slide*N*N+row*N+col];
            continue;
          }
        
          // get current temperature at (x,y,z)
          value_t tc = A[row*N+col+slide*N*N];

          // get temperatures left/right and up/down

          value_t tl = (col != 0    ) ? A[slide*N*N + row*N + (col-1)] : tc;
          value_t tr = (col != N - 1) ? A[slide*N*N + row*N + (col+1)] : tc;
          value_t tu = (row != 0    ) ? A[slide*N*N + (row-1)*N + col] : tc;
          value_t td = (row != N - 1) ? A[slide*N*N + (row+1)*N + col] : tc;
          value_t tb = (slide != 0    ) ? A[(slide-1)*N*N + row*N + col] : tc;
          value_t ta = (slide != N - 1) ? A[(slide+1)*N*N + row*N + col] : tc;

          // update temperature at current point

          value_t temp  = tc + 1.0/7 * (tl + tr + tu + td + tb + ta + (-6*tc));
          B[row*N+col+slide*N*N] = temp;

        }
      }
    }
    // swap matrices (just pointers, not content)
    Vector H = A;
    A = B;
    B = H;
    
    // show intermediate step
    if (!(t % 1000) && 0) {
      printf("Step t=%d:\n", t);
      printTemperature(A, N, N, N);
      printf("\n");
    }
  }
  

  subsize[0] = partsize_row;
  subsize[1] = partsize_col;
  subsize[2] = partsize_slide;
  
  substart[0] = start_row;
  substart[1] = start_col;
  substart[2] = start_slide;
  int sendcount[size];
  int displacement[size];
  for (int i=0; i < size; ++i) {
    sendcount[i]=1;
  }
  int colcount = 0;
  int rowcount = 0;
  for (int i=0; i < size; ++i) {

    displacement[i]=colcount + rowcount;
    colcount += partsize_col;
    if ((i+1) % dim[1] == 0){
      colcount = 0;
      rowcount += partsize_row * N;
    }
  }
  
  /* printf("subsize %d %d start %d %d\n",subsize[0], subsize[1], substart[0], substart[1] ); */
  /* MPI_Type_create_subarray(DIMENSIONS, arraysize, subsize, substart, MPI_ORDER_C, MPI_DOUBLE, &gather); */
  /* MPI_Type_commit(&gather); */
  /* MPI_Type_create_resized(gather, 0, 1*sizeof(double), &gather2); */
  /* MPI_Type_commit(&gather2); */

  /* MPI_Gatherv(A, 1, gather, A, sendcount, displacement, gather2, 0, MPI_COMM_WORLD); */
  /* MPI_Type_free(&gather); */
  /* MPI_Type_free(&gather2); */

// ---------- check ----------

  

  int success = 1;
  printf("Final:\n");
  printTemperature(A, N, N, N);
  printf("\n");
  
  for (int z = 0; z < N; z++) {
    for (int y = 0; y < N; y++) {
      for (int x = 0; x < N; x++) {
        value_t temp = A[z*N*N + y*N + x];
        if (273 <= temp && temp <= 273 + 60)
          continue;
        success = 0;
        printf("failure on cell x=%d y=%d z=%d with value %f\n", x, y, z, temp);
        break;
      }
    }
  }
  printf("Verification: %s\n", (success) ? "OK" : "FAILED");
  



  // ---------- cleanup ----------

  releaseVector(A);
  releaseVector(B);


  MPI_Finalize();

  // done
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Vector createVector(int N, int M, int Oz) {
  // create data and index vector
  return malloc(sizeof(value_t) * N * M * Oz);
}

void releaseVector(Vector m) { free(m); }

// taken from old parallel course https://git.uibk.ac.at/csat2062/parallel_local
void printTemperature(Vector m, int N, int M, int Oz) {
  const char *colors = " .-:=+*#%@";
  const int numColors = 10;

  // boundaries for temperature (for simplicity hard-coded)
  const value_t max = 273 + 30;
  const value_t min = 273 + 0;

  // set the 'render' resolution
  int H = 30;
  int W = 50;
  int D = 50;

  // step size in each dimension
  int sH = N / H;
  int sW = M / W;
  int sD = Oz / D;

  // upper wall
  for (int i = 0; i < W + 2; i++) {
    printf("X");
  }
  printf("\n");

  // room
  for (int i = 0; i < H; i++) {
    // left wall
    printf("X");
    // actual room
    for (int j = 0; j < W; j++) {

      // get max temperature in this tile
      value_t max_t = 0;
      for (int x = sH * i; x < sH * i + sH; x++) {
        for (int y = sW * j; y < sW * j + sW; y++) {
          //add values over z coordinate
          for (int z = 0; z < Oz; z++) {
            max_t = (max_t < m[x * N + y + N*M*z]) ? m[x * N + y+ N*M*z] : max_t;
          }
        }
      }
      // avarage over z hight
      value_t temp = max_t;

      // pick the 'color'
      int c = ((temp - min) / (max - min)) * numColors;
      c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);

      // print the average temperature
      printf("%c", colors[c]);
    }
    // right wall
    printf("X\n");
  }

  // lower wall
  for (int i = 0; i < W + 2; i++) {
    printf("X");
  }
  printf("\n");
}

void neighbourMPI(int dimension, int* arraysize, int* subsize, int* substart, int* neighbour, MPI_Comm cartesian, double* buffer, enum mode mode){
  int neighbourrank;
  MPI_Datatype subarray;
  MPI_Request req;
  MPI_Cart_rank(cartesian, neighbour, &neighbourrank);
  //printf("mypos: %d %d %d creating subarray @ %d %d %d with size %d %d %d mode %d for/from neighbour %d %d %d with rank %d \n",coord[0], coord[1], coord[2], substart[0], substart[1], substart[2], subsize[0], subsize[1], subsize[2], mode, neighbour[0], neighbour[1], neighbour[2], neighbourrank );
  //only create 2d subarrays
  MPI_Type_create_subarray(dimension-1, arraysize, subsize, substart, MPI_ORDER_C, MPI_DOUBLE, &subarray);
  MPI_Type_commit(&subarray);
  //send one by one
  if (mode == SEND){
    for (int i = 0; i < subsize[2]; ++i) {
      MPI_Send(buffer+(substart[2]+i)*arraysize[0]*arraysize[1], 1, subarray, neighbourrank, 0, MPI_COMM_WORLD);
    }
    //printf("mypos: %d %d %d send value %f\n", coord[0], coord[1], coord[2], buffer[substart[0]*arraysize[0] + substart[1]+ substart[2]*arraysize[0]*arraysize[0]]);
  } else {
    for (int i=0; i < subsize[2]; ++i) {
      MPI_Recv(buffer+(substart[2]+i)*arraysize[0]*arraysize[1], 1, subarray, neighbourrank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    //printf("mypos: %d %d %d value after receive %f for position %d %d %d\n", coord[0], coord[1], coord[2], buffer[substart[0]*arraysize[0] + substart[1]+ substart[2]*arraysize[0]*arraysize[0]], substart[0], substart[1], substart[2]);
  }
  MPI_Type_free(&subarray);

}
