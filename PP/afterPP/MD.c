/*  Simple molecular dynamics code.   */

#include <stdio.h>
#include <math.h>
#include "coord.h"

void evolve( int count ){
int  step, i, j, k, l, collided;
double Size;

/*  Loop over timesteps  */
    for( step = 1 ; step <= count ; step++ ){

        printf( "timestep %d\n", step );
        printf( "collisions %d\n", collisions );

/* set the viscosity and the wind term in the force calculation */
        for( i = 0 ; i < Nbody ; i++ ){

            f[0][i] = -vis[i] * velo[0][i] -vis[i] * wind[0];
            f[1][i] = -vis[i] * velo[1][i] -vis[i] * wind[1];
            f[2][i] = -vis[i] * velo[2][i] -vis[i] * wind[2];      
        }

/* calculate central force */
        for( i = 0 ; i < Nbody ; i++ ){

 	  for( l = 0 ; l < Ndim ; l++ ){

                f[l][i] = f[l][i] - force( G * mass[i] * M_central, pos[l][i], sqrt( (pos[0][i] * pos[0][i]) + (pos[1][i] * pos[1][i]) + (pos[2][i] * pos[2][i]) ));
	  }
	}

/*  add pairwise forces. */

        k = 0;

        for( i = 0 ; i < Nbody ; i++ ){
          
          for( j = i + 1 ; j < Nbody ; j++ ){

            Size = radius[i] + radius[j];

            collided=0;

            delta_pos[0][k] = pos[0][i] - pos[0][j];

            delta_pos[1][k] = pos[1][i] - pos[1][j];

            delta_pos[2][k] = pos[2][i] - pos[2][j];
   
            double r2 = (delta_pos[0][k] * delta_pos[0][k]) + (delta_pos[1][k] * delta_pos[1][k]) + (delta_pos[2][k] * delta_pos[2][k]) ;
            
            double r3 = sqrt(r2);
   
/*  flip force if close in */

              if( r3 >= Size ){

                for( l = 0 ; l < Ndim ; l++ ){

                    double temp01 = force( G * mass[i] * mass[j], delta_pos[l][k], r3);

                    f[l][i] = f[l][i] - temp01;

                    f[l][j] = f[l][j] + temp01;

                 }
              }else{

                for( l = 0 ; l < Ndim ; l++ ){

                    double temp02 = force( G * mass[i] * mass[j], delta_pos[l][k], r3);
 
                    f[l][i] = f[l][i] + temp02;

                    f[l][j] = f[l][j] - temp02;

	   	collided=1;

              }
         
            }

	    if( collided == 1 ){

	      collisions++;

	    }

            k = k + 1;

          }

        }
/* update the velocity and the positions */
        for(j=0;j<Ndim;j++){ 

            for( i = 0 ; i < Nbody ; i++ ){

                pos[j][i] = pos[j][i] + dt * velo[j][i];

                velo[j][i] = velo[j][i] + dt * (f[j][i]/mass[i]); 
            }
        }

      }

}
