#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

const int DEBUGPRINT = 0;
const float G = 1.0f;
const float EPS = 0.0001f;

struct body {
  float pos[2];
  float vel[2];
  float mass;
};

void printbody(struct body *a) {
  printf("pos=(%f, %f) vel=(%f, %f) mass=%f\n", a->pos[0], a->pos[1], a->vel[0],
         a->vel[1], a->mass);
}

// random float between low and high
float randomfloat(unsigned int *seed, float low, float high) {
  return (float)rand_r(seed) / (float)RAND_MAX * (high - low) + low;
}

void randomize_body(struct body *a, unsigned int *seed) {
  a->pos[0] = randomfloat(seed, -1000.0f, 1000.0f);
  a->pos[1] = randomfloat(seed, -1000.0f, 1000.0f);
  a->vel[0] = randomfloat(seed, -1.0f, 1.0f);
  a->vel[1] = randomfloat(seed, -1.0f, 1.0f);
  a->mass = randomfloat(seed, 1.0f, 1000.0f);
}

void calculate_force_float(struct body *a, struct body *b, float *returnvalue) {
  float distancex = b->pos[0] - a->pos[0];
  float signx = (distancex < 0.0f) ? -1.0f : 1.0f;
  float distancey = b->pos[1] - a->pos[0];
  float signy = (distancey < 0.0f) ? -1.0f : 1.0f;
  float distancesquared = distancex * distancex + distancey * distancey;
  float forcediag =
      G * a->mass * b->mass / powf((distancesquared / 4.0f + EPS), 3.0f / 2.0f);
  returnvalue[0] =
      forcediag * signx * (distancex * distancex) / distancesquared;
  returnvalue[1] =
      forcediag * signy * (distancey * distancey) / distancesquared;
}

void updateforce_float(struct body *a, float forcex, float forcey, int dt) {
  a->vel[0] += forcex / a->mass * dt;
  a->vel[1] += forcey / a->mass * dt;
}

void updateforce_float2(struct body *a, float *force, int dt) {
  a->vel[0] += force[0] / a->mass * dt;
  a->vel[1] += force[1] / a->mass * dt;
}

void movebody(struct body *a, int dt) {
  a->pos[0] += a->vel[0] * dt;
  a->pos[1] += a->vel[1] * dt;
}

int main(int argc, char *argv[]) {
  int N = 10000;
  int T = 100;
  if (argc == 3) {
    N = atoi(argv[1]); // number of particles
    T = atoi(argv[2]); // number of timesteps
  }
  srand(1);
  struct body *bodyarray = malloc(sizeof(struct body) * N);
#pragma omp parallel
  {
    printf("number of threads =%d\n", omp_get_num_threads());
    unsigned int seed = omp_get_thread_num();
#pragma omp for
    for (int i = 0; i < N; ++i) {
      randomize_body(&bodyarray[i], &seed);
      if (DEBUGPRINT)
        printbody(&bodyarray[i]);
    }
  }

  float forceupdate[2];
  for (int t = 0; t < T; ++t) {
#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
      for (int j = i + 1; j < N; ++j) {
        calculate_force_float(&bodyarray[i], &bodyarray[j], forceupdate);
        if (DEBUGPRINT)
          printf("force=(%f, %f)\n", forceupdate[0], forceupdate[1]);
        updateforce_float2(&bodyarray[i], forceupdate, 1);
        forceupdate[0] *= -1.0f;
        forceupdate[1] *= -1.0f;
        updateforce_float2(&bodyarray[j], forceupdate, 1);  //race condition possible
      }
    }
#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
      movebody(&bodyarray[i], 1);
    }
  }
  if (DEBUGPRINT) {
    printf("****************\n");
    for (int i = 0; i < N; ++i) {
      printbody(&bodyarray[i]);
    }
  }

  free(bodyarray);
  return 0;
}
