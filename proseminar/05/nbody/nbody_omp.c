#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define NDEBUG

const float G = 1.0f;
const float EPS = 0.0001f;

struct f_touple {
  float x;
  float y;
};

struct body {
  struct f_touple pos;
  struct f_touple vel;
  float mass;
};

void printbody(struct body *a) {
  printf("pos=(%f, %f) vel=(%f, %f) mass=%f\n", a->pos.x, a->pos.y, a->vel.x,
         a->vel.y, a->mass);
}

// random float between low and high
float randomfloat(unsigned int *seed, float low, float high) {
  return (float)rand_r(seed) / (float)RAND_MAX * (fabsf(low - high)) + low;
}

void randomize_body(struct body *a, unsigned int *seed) {
  a->pos.x = randomfloat(seed, -1000.0, 1000.0);
  a->pos.y = randomfloat(seed, -1000.0, 1000.0);
  a->vel.x = randomfloat(seed, -1.0, 1.0);
  a->vel.y = randomfloat(seed, -1.0, 1.0);
  a->mass = randomfloat(seed, 0.0, 1000.0);
}

struct f_touple calculate_force(struct body *a, struct body *b) {
  float distancex = b->pos.x - a->pos.x;
  float distancey = b->pos.y - a->pos.y;
  float distance = sqrt(distancex*distancex+distancey*distancey);
  struct f_touple force;
  float forcediag = G * a->mass * b->mass /
            powf(((distance / 2) * (distance / 2) + EPS), 3.0 / 2.0);
  force.x = forcediag * distancex / distance;
  force.y = forcediag * distancey / distance;
  return force;
}

void updateforce(struct body *a, struct f_touple force, int dt) {
  a->vel.x += force.x / a->mass * dt;
  a->vel.y += force.y / a->mass * dt;
}

void movebody(struct body *a, int dt) {
  a->pos.x += a->vel.x * dt;
  a->pos.y += a->vel.y * dt;
}

int main(int argc, char *argv[]) {
  int N = 10000;
  int T = 1000;
  if (argc == 3) {
    N = atoi(argv[1]);
    T = atoi(argv[2]);
  }
  struct body *test = malloc(sizeof(struct body) * N);
  #pragma omp parallel for
  for (int i = 0; i < N; ++i) {
    unsigned int seed = omp_get_thread_num()+i;
    randomize_body(&test[i], &seed);
#ifndef NDEBUG
    printbody(&test[i]);
#endif
  }

  for (int t = 0; t < T; ++t) {
#pragma omp parallel for collapse(1)
    for (int i = 0; i < N; ++i) {
      for (int j = i + 1; j < N; ++j) {
        #pragma omp critical
        {
        struct f_touple force = calculate_force(&test[i], &test[j]);
        //printf("force=(%f, %f)\n", force.x, force.y);
        updateforce(&test[i], force, 1);
        force.x *= -1;
        force.y *= -1;
        updateforce(&test[j], force, 1);
        }
      }
    }
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
      movebody(&test[i], 1);
    }
  }
#ifndef NDEBUG
  printf("****************\n");
  for (int i = 0; i < N; ++i) {
    printbody(&test[i]);
  }
#endif
  free(test);
  return 0;
}
