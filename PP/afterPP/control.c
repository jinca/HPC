/*
 *
 * Control program for the MD update
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#define DECL
#include "coord.h"

double second(void); 
int main(int argc, char *argv[]){
  int i,j;
  FILE *in, *out;
  double tstart,tstop;
  double start,stop;
  int Nstep=100;
  char name[80];

  if( argc > 1 ){
    Nstep = atoi( argv[1] );
  }
  wind[Xcoord] = 0.9;
  wind[Ycoord] = 0.4;
  wind[Zcoord] = 0.0;
/*  
    r = malloc(Nbody*sizeof(double)); 
    delta_r = malloc(Npair*sizeof(double)); 
    mass = malloc(Nbody*sizeof(double));
    radius = malloc(Nbody*sizeof(double));
    vis = malloc(Nbody*sizeof(double));
    for(i=0;i<Ndim;i++){
        f[i] = malloc(Nbody*sizeof(double));
        pos[i] = malloc(Nbody*sizeof(double));
        velo[i] = malloc(Nbody*sizeof(double));
        delta_pos[i] = malloc(Npair*sizeof(double));
    }
*/
/* read the initial data from a file */

  collisions = 0;
  in = fopen( "input.dat", "r" );

  if( ! in ) 
  {
      perror( "input.dat" );
      exit(1);
  }

  for(i=0;i<Nbody;i++){
    fscanf(in,"%16le%16le%16le%16le%16le%16le%16le%16le%16le\n",
      mass+i,radius+i,vis+i,
      &pos[Xcoord][i], &pos[Ycoord][i], &pos[Zcoord][i],
      &velo[Xcoord][i], &velo[Ycoord][i], &velo[Zcoord][i]);
  }
  fclose(in);

/* Run Nstep timesteps and time how long it took  */
 
   tstart = second();
   for( j = 1 ; j <= Nsave ; j++ ){
      start = second();
      evolve( Nstep ); 
      stop = second();
      printf( "%d timesteps took %f seconds\n", Nstep, stop-start);
      printf( "collisions %d\n", collisions);

/* write final result to a file */
      sprintf(name,"output.dat%03d",j*Nstep);
      out = fopen(name,"w");

      if( ! out ){
	perror(name);
	exit(1);
      }

      for(i=0;i<Nbody;i++){
	fprintf(out,"%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E\n",
		mass[i],radius[i],vis[i],
		pos[Xcoord][i], pos[Ycoord][i], pos[Zcoord][i],
		velo[Xcoord][i], velo[Ycoord][i], velo[Zcoord][i]);
      }
      fclose(out);
  }
  tstop=second();
  printf("%d timesteps took %f seconds\n",Nsave*Nstep,tstop-tstart);

}

double second()
{
        struct timeval tp;
        struct timezone tzp;
        int i;

        i = gettimeofday(&tp,&tzp);
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}
