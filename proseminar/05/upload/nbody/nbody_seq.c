#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const int DEBUGPRINT = 0;
const float G = 1.0f;
const float EPS = 0.0001f;

struct f_touple
{
  float x;
  float y;
};

struct body
{
  float pos[2];
  float vel[2];
  float mass;
};

void printbody(struct body *a)
{
  printf("pos=(%f, %f) vel=(%f, %f) mass=%f\n", a->pos[0], a->pos[1], a->vel[0],
         a->vel[1], a->mass);
}

// random float between low and high
float randomfloat(unsigned int *seed, float low, float high)
{
  return (float)rand_r(seed) / (float)RAND_MAX * (high - low) + low;
}

void randomize_body(struct body *a, unsigned int *seed)
{
  a->pos[0] = randomfloat(seed, -1000.0, 1000.0);
  a->pos[1] = randomfloat(seed, -1000.0, 1000.0);
  a->vel[0] = randomfloat(seed, -1.0, 1.0);
  a->vel[1] = randomfloat(seed, -1.0, 1.0);
  a->mass = randomfloat(seed, 0.0, 1000.0);
}

struct f_touple calculate_force(struct body *a, struct body *b)
{
  float distancex = b->pos[0] - a->pos[0];
  float signx = (distancex < 0.0) ? -1.0 : 1.0;
  float distancey = b->pos[1] - a->pos[1];
  float signy = (distancey < 0.0) ? -1.0 : 1.0;
  float distancesquared = distancex * distancex + distancey * distancey;
  struct f_touple force;
  float forcediag = G * a->mass * b->mass /
                    powf((distancesquared / 4.0 + EPS), 3.0 / 2.0);
  force.x = forcediag * signx * (distancex * distancex) / distancesquared;
  force.y = forcediag * signy * (distancey * distancey) / distancesquared;
  return force;
}

void calculate_force_float(struct body *a, struct body *b, float *returnvalue)
{
  float distancex = b->pos[0] - a->pos[0];
  float signx = (distancex < 0.0) ? -1.0 : 1.0;
  float distancey = b->pos[1] - a->pos[0];
  float signy = (distancey < 0.0) ? -1.0 : 1.0;
  float distancesquared = distancex * distancex + distancey * distancey;
  float forcediag = G * a->mass * b->mass /
                    powf((distancesquared / 4.0 + EPS), 3.0 / 2.0);
  returnvalue[0] = forcediag * signx * (distancex * distancex) / distancesquared;
  returnvalue[1] = forcediag * signy * (distancey * distancey) / distancesquared;
}

void updateforce(struct body *a, struct f_touple force, int dt)
{
  a->vel[0] += force.x / a->mass * dt;
  a->vel[1] += force.y / a->mass * dt;
}

void updateforce_float(struct body *a, float forcex, float forcey, int dt)
{
  a->vel[0] += forcex / a->mass * dt;
  a->vel[1] += forcey / a->mass * dt;
}

void movebody(struct body *a, int dt)
{
  a->pos[0] += a->vel[0] * dt;
  a->pos[1] += a->vel[1] * dt;
}

int main(int argc, char *argv[])
{
  int N = 10000;
  int T = 100;
  if (argc == 3)
  {
    N = atoi(argv[1]);
    T = atoi(argv[2]);
  }
  srand(1);
  struct body *test = malloc(sizeof(struct body) * N);
  unsigned int seed = 1;
  for (int i = 0; i < N; ++i)
  {
    seed = i;
    randomize_body(&test[i], &seed);
    if (DEBUGPRINT)
      printbody(&test[i]);
  }
  //struct f_touple force;
  float forceupdate[2];
  for (int t = 0; t < T; ++t)
  {

    for (int i = 0; i < N; ++i)
    {
      for (int j = i + 1; j < N; ++j)
      {
        //force = calculate_force(&test[i], &test[j]);
        calculate_force_float(&test[i], &test[j], forceupdate);
        if (DEBUGPRINT)
          printf("force=(%f, %f)\n", forceupdate[0], forceupdate[1]);
        //updateforce(&test[i], force, 1);
        updateforce_float(&test[i], forceupdate[0], forceupdate[1], 1);
        //force.x *= -1.0;
        forceupdate[0] *= -1.0;
        //force.y *= -1.0;
        forceupdate[1] *= -1.0;
        //updateforce(&test[j], force, 1);
        updateforce_float(&test[j], forceupdate[0], forceupdate[1], 1);
      }
    }
    for (int i = 0; i < N; ++i)
    {
      movebody(&test[i], 1);
    }
  }
  if (DEBUGPRINT)
  {
    printf("****************\n");
    for (int i = 0; i < N; ++i)
    {
      printbody(&test[i]);
    }
  }

  free(test);
  return 0;
}
