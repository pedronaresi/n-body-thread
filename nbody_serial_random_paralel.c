#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>

#define MODULUS    2147483647
#define MULTIPLIER 48271
#define DEFAULT    123456789

#define MAX_THREADS 2 //1, 2, 4, 8

typedef struct {
  double x, y, z;
  double mass;
} Particle;
typedef struct {
  double xold, yold, zold;
  double fx, fy, fz;
} ParticleV;

double Random(int thid);
void* InitParticles(void *);
void* ComputeForces(void *);
void* ComputeNewPos(void *);

Particle  * particles;   /* Particulas */
ParticleV * pv;          /* Velocidade da Particula */
int         npart;
double    * sim_t;       /* Tempo de Simulacao */

long *seed;

struct timeval inicio, final2;
int tmili;

double * max_f;

int main() {
  int th;
  //scanf("%d",&npart);
  npart = 25000;
  /* Allocate memory for particles */
  particles = (Particle *) malloc(sizeof(Particle)*npart);
  pv = (ParticleV *) malloc(sizeof(ParticleV)*npart);
  sim_t = (double *) malloc(sizeof(double)*MAX_THREADS);
  max_f = (double *) malloc(sizeof(double)*MAX_THREADS);
  seed = (long *) malloc(sizeof(long)*MAX_THREADS);
  /* Generate the initial values */
  
  //Creating array of Threads
  pthread_t t[MAX_THREADS];

  for (int i = 0; i < MAX_THREADS; i++) {
    sim_t[i] = 0.0;
    seed[i]=DEFAULT;
  }

  gettimeofday(&inicio, NULL);
  
  for(th=0; th<MAX_THREADS; th++) {
    pthread_create(&t[th], NULL, InitParticles, (void *) th);
  }

  for(th=0; th<MAX_THREADS; th++) {
    pthread_join(t[th],NULL);
  }
  
  for(th=0; th<MAX_THREADS; th++) {
    pthread_create(&t[th], NULL, ComputeForces, (void *) th);
  }

  for(th=0; th<MAX_THREADS; th++) {
    pthread_join(t[th],NULL);
  }

  for(th=1; th<MAX_THREADS; th++) {
    if (max_f[0]<max_f[th]) max_f[0]=max_f[th];
  }

  for(th=0; th<MAX_THREADS; th++) {
    pthread_create(&t[th], NULL, ComputeNewPos, (void *) th);
  }

  for(th=0; th<MAX_THREADS; th++) {
    pthread_join(t[th],NULL);
  }
  
  gettimeofday(&final2, NULL);
  tmili = (int) (1000 * (final2.tv_sec - inicio.tv_sec) +
		 (final2.tv_usec - inicio.tv_usec) / 1000);

  //for (i=0; i<npart; i++)
  //fprintf(stdout,"%.5lf %.5lf %.5lf\n", particles[i].x, particles[i].y, particles[i].z);

  printf("%g\n", max_f[0]);
  printf("%g\n", sim_t[0]);

  printf("%d\n", tmili);

  return 0;
}

double Random(int thid) {
  const long Q = MODULUS / MULTIPLIER;
  const long R = MODULUS % MULTIPLIER;
  long t;

  t = MULTIPLIER * (seed[thid] % Q) - R * (seed[thid] / Q);
  if (t > 0)
    seed[thid] = t;
  else
    seed[thid] = t + MODULUS;
  return ((double) seed[thid] / MODULUS);
}

void *InitParticles(void *tid) {
  int thid;
  thid = (long) tid;
  
  printf("Init_Particles thread: %d\n", thid);

  for (int i=thid; i<npart; i+=MAX_THREADS) {
    particles[i].x	  = Random(thid);
    particles[i].y	  = Random(thid);
    particles[i].z	  = Random(thid);
    particles[i].mass = 1.0;
    pv[i].xold	  = particles[i].x;
    pv[i].yold	  = particles[i].y;
    pv[i].zold	  = particles[i].z;
    pv[i].fx	  = 0;
    pv[i].fy	  = 0;
    pv[i].fz	  = 0;
  }

  pthread_exit(0);
}

void *ComputeForces(void *tid) {
  double max_in_f;
  int i, j;
  max_in_f = 0.0;

  int thid;

  thid = (long) tid;

  printf("  Compute_Forces thread: %d\n", thid);

  double xi, yi, rx, ry, mj, r, fx, fy, rmin;

  for (i=thid; i<npart; i+=MAX_THREADS) {
    rmin = 100.0;
    xi   = particles[i].x;
    yi   = particles[i].y;
    fx   = 0.0;
    fy   = 0.0;
    for (j=thid; j<npart; j+=MAX_THREADS) {
      rx = xi - particles[j].x;
      ry = yi - particles[j].y;
      mj = particles[j].mass;
      r  = rx * rx + ry * ry;
      /* ignore overlap and same particle */
      if (r == 0.0) continue;
      if (r < rmin) rmin = r;
      r  = r * sqrt(r);
      fx -= mj * rx / r;
      fy -= mj * ry / r;
    }
    pv[i].fx += fx;
    pv[i].fy += fy;
    fx = sqrt(fx*fx + fy*fy)/rmin;
    if (fx > max_in_f) max_in_f = fx;
  }
  max_f[thid] = max_in_f;

  pthread_exit(0);
}

void *ComputeNewPos(void *tid) {
  int i;
  double a0, a1, a2;
  static double dt_old = 0.001, dt = 0.001;
  double dt_new;

  int thid;
  thid = (long) tid;

  printf("    Compute_New_Pos thread: %d\n", thid);

  a0	 = 2.0 / (dt * (dt + dt_old));
  a2	 = 2.0 / (dt_old * (dt + dt_old));
  a1	 = -(a0 + a2);

  for (i=thid; i<npart; i+=MAX_THREADS) {
    double xi, yi;
    xi	           = particles[i].x;
    yi	           = particles[i].y;
    particles[i].x = (pv[i].fx - a1 * xi - a2 * pv[i].xold) / a0;
    particles[i].y = (pv[i].fy - a1 * yi - a2 * pv[i].yold) / a0;
    pv[i].xold     = xi;
    pv[i].yold     = yi;
    pv[i].fx       = 0;
    pv[i].fy       = 0;
  }
  dt_new = 1.0/sqrt(max_f[0]);
  /* Set a minimum: */
  if (dt_new < 1.0e-6) dt_new = 1.0e-6;
  /* Modify time step */
  if (dt_new < dt) {
    dt_old = dt;
    dt     = dt_new;
  }
  else if (dt_new > 4.0 * dt) {
    dt_old = dt;
    dt    *= 2.0;
  }
  sim_t[thid] = dt_old;
  pthread_exit(0);
}
