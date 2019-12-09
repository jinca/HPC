/* MPP Coursework by B129706 
 * Parallel code from the serial version: imagenew.c  
 * Use for this exercise 
   -> edgenew768x768.pgm 
   -> pgmio.c
   -> pgmio.h
   -> MAXITER 1000
   It will produce imagenewp768x768_1000_test.pgm          */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include "pgmio.h"
#define M 768
#define N 768
#define SHIFT_ROW 0
#define SHIFT_COL 1
#define DISP 1
#define MAXITER 1000 
#define PRINTFREQ  200

double boundaryval( int i, int m );

int main( int argc, char *argv[] )
{
   /* variables to use an MPI environment */
   int         size;
   int         rank;
   int         source; 
   int         tag;
   MPI_Status  status; 
   MPI_Request request
;
   /* variables to create a topology */
   int         ndims = 2;
   int         reorder;
   int         cart_rank;
   int         nbr_down, nbr_up, nbr_left, nbr_right;
   int         dims[ndims]; dims[0]=0; dims[0]=0;
   int         coord[ndims];
   int         period[ndims];
   MPI_Comm    comm2D;

   /* variables to generate the image */
   double      buf[M][N];
   double      val;
   int         i,j,iter;
   char        *filename;

   /* variables to meassure time of the program */
   double      startim;
   double      endtim;
   double      total;

   /* start up MPI environment: MPI_COMM_WORLD */
   MPI_Init( &argc, &argv );
   MPI_Comm_size( MPI_COMM_WORLD, &size );
   MPI_Comm_rank( MPI_COMM_WORLD, &rank );

   /* create cartesian topology for processes */
   MPI_Dims_create( size, ndims, dims );
   if( rank == 0 )
   {
       printf( "Currently working [%d]/[%d]", rank, size );
       printf( "Dimensions = [%d x %d] \n", dims[0], dims[1] );
       printf( "Processing %d x %d image\n", M, N );
       printf( "Number of iterations = %d\n", MAXITER );
       filename = "edgenew256x192.pgm";
       printf( "\n Reading <%s>\n", filename );
       pgmread( filename, buf, M, N );
       printf( "\n" );
   }

   /* measure time of this parallel code */
   startim = MPI_Wtime();
   MPI_Barrier(MPI_COMM_WORLD);

   /* broadcast the image in MPI_COMM_WORLD */ 
   MPI_Bcast( &buf, M*N, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  
   /****** Test 1 *****/ 
   if ( rank == 1 )
         pgmwrite( "afterBCast.pgm", buf, M, N );

   /* create cartesian mapping */
   period[0] = 0;
   period[1] = 1;
   reorder = 0;
   MPI_Cart_create( MPI_COMM_WORLD, ndims, dims, period, reorder, &comm2D );

   MPI_Cart_coords( comm2D, rank, ndims, coord );
   MPI_Cart_rank( comm2D, coord, &cart_rank );
 
   MPI_Cart_shift( comm2D, SHIFT_ROW, DISP, &nbr_left, &nbr_right );
   MPI_Cart_shift( comm2D, SHIFT_COL, DISP, &nbr_down, &nbr_up );

   /****** Test 2 *****/ 
   printf( "Rank[%d] Coord(%d,%d): SHIFT_DIM[%d],Shift=%d:", 
          rank, coord[0], coord[1], SHIFT_ROW, DISP ); 
   printf( "\n nbr_left[%d] P[%d] nbr_right[%d] \n", 
          nbr_left, rank, nbr_right );
   printf( "Rank[%d] Coord(%d,%d): SHIFT_DIM[%d],Shift=%d:", 
          rank, coord[0], coord[1], SHIFT_COL, DISP );
   printf( "\n nbr_down[%d] P[%d] nbr_up[%d] \n", 
          nbr_down, rank, nbr_up );

   /* declare image chunks variables */
   int Mp = M/dims[0];
   int Np = N/dims[1];
   int x_coord, y_coord;
   double old[Mp+2][Np+2], new[Mp+2][Np+2], edge[Mp+2][Np+2];

   /* get the values of the coordinates */
   x_coord = coord[0]*Mp;
   y_coord = coord[1]*Np;

   /* variable to be used in the halos exchange */
   MPI_Datatype row;
   MPI_Type_vector( Mp,1,Np+2, MPI_DOUBLE, &row );
   MPI_Type_commit( &row );
   
   /****** Test 3 *****/ 
   printf( "Rank [%d] has Mp: %d Np: %d x_coord: %d, y_coord: %d",rank,Mp,Np,x_coord,y_coord );
  
   /* generate the image from the edges */
   for( i = 1 ; i < Mp+1 ; i++ )
   {
       for( j = 1 ; j < Np+1 ; j++ )
       {
           edge[i][j] = buf[x_coord + i-1][y_coord + j-1];
       }
   }
 
   for( i = 0 ; i < Mp+2 ; i++ )
   {
       for ( j=0 ; j < Np+ 2; j++ )
       {
          old[i][j] = 255.0;
       }
   }

  /* set fixed boundary conditions on the left and right sides */

   for( j = 1 ; j < N ; j++ )
   {
      /* compute sawtooth value */

      val = boundaryval( j, N );

      old[0][j] = (int)(255.0*(1.0-val));
      old[M][j] = (int)(255.0*val);
   }

   for( iter = 1 ; iter <= MAXITER; iter++ )
   {
      if( iter % PRINTFREQ == 0 )
      {
          printf("Iteration %d\n", iter);
      }

      /* Implement periodic boundary conditions on bottom and top sides */

      for( i = 1 ; i < Mp+1 ; i++ )
      {
          old[i][0]   = old[i][Np+1];
          old[i][Np+1] = old[i][1];
      }

      for( i = 1 ; i < Mp+1 ; i++ )
      {
          for ( j = 1; j < Np+1 ; j++)
          {
              new[i][j]=0.25*(old[i-1][j]+old[i+1][j]+old[i][j-1]+old[i][j+1] - edge[i][j]);
          }
      }

      for( i = 1 ; i < Mp+1 ; i++ )
      {
          for ( j = 1 ; j< Np+1 ; j++ )
          {
              old[i][j] = new[i][j];
          }
      }

       /* send rows to up neighbor */
       MPI_Isend( &old[1][Np], 1, row, nbr_up, 11, MPI_COMM_WORLD,&request );
       MPI_Recv( &old[1][Np+1],1, row, nbr_down, 11, MPI_COMM_WORLD, &status );
       MPI_Wait( &request, MPI_STATUS_IGNORE );

       /* send rows to down neighbor */
       MPI_Isend( &old[1][1], 1, row, nbr_down, 22, MPI_COMM_WORLD,&request );
       MPI_Recv( &old[1][0],1,row, nbr_up, 22, MPI_COMM_WORLD, &status );
       MPI_Wait( &request, MPI_STATUS_IGNORE );

       /* send column to right neighbor */
       MPI_Isend( &old[Mp][1], Np, MPI_DOUBLE, nbr_right, 33, MPI_COMM_WORLD, &request );
       MPI_Recv( &old[0][1], Np, MPI_DOUBLE, nbr_left, 33, MPI_COMM_WORLD, &status );
       MPI_Wait( &request, MPI_STATUS_IGNORE );

       /* send column to left neighbor */
       MPI_Isend( &old[1][1], Np, MPI_DOUBLE, nbr_left, 44, MPI_COMM_WORLD, &request );
       MPI_Recv( &old[Mp+1][1], Np, MPI_DOUBLE, nbr_right,44, MPI_COMM_WORLD, &status );
       MPI_Wait( &request, MPI_STATUS_IGNORE );
  
    }

    printf("\nFinished %d iterations\n", iter-1);

    /* send the chunks of the image from different processors to rank 0 */
    if ( rank == 0 )
    {
       /****** Test 4 *****/ 
       printf( "\nRank %d has coords: (%d,%d)\n", rank, coord[0], coord[1] );
      
       /* map the first chunk of the image */ 
       for( i = 1; i < Mp+1; i++ )
       {
          for( j = 1; j < Np+1; j++) 
          { 
             buf[x_coord+i-1][y_coord+j-1] = old[i][j];
          }
       }  
      
       /* receive the chunks from processors different of zero */ 
       for ( source = 1; source < size ;  source ++)
       {
          tag = 0;
          MPI_Recv( old,(Mp+2)*(Np+2), MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status );
          MPI_Recv( coord, 2, MPI_INT, source, tag, MPI_COMM_WORLD, &status );

          /****** Test 5 *****/ 
          printf( "Rank 0 receives from rank %d: coords = %d %d\n", source, coord[0], coord[1] );

          /* map the chunks of the rest processors */
          x_coord = coord[0]*Mp;
          y_coord = coord[1]*Np;
        
          for( i = 1; i < Mp+1; i++ )
          {
              for( j = 1; j < Np+1; j++) 
              { 
                  buf[x_coord+i-1][y_coord+j-1] = old[i][j];
              }
          }  
      } 
   }
   else
   {
      tag = 0;
      /****** Test 6 *****/ 
      printf( "Rank %d sends to rank 0\n", rank );

      /* send the chunks of the image to rank zero */
      MPI_Send( old, (Mp+2)*(Np+2), MPI_DOUBLE, 0, tag, MPI_COMM_WORLD );
      MPI_Send( coord, 2, MPI_INT, 0, tag, MPI_COMM_WORLD );
   }

   /* calculate elapse time */
   MPI_Barrier( MPI_COMM_WORLD );
   endtim = MPI_Wtime() - startim;
  
   if ( rank == 0 )
   {
      /* print the chunks of the image into buf */
      filename = "imagenew256x192_test.pgm";
      printf( "\nWriting <%s>\n", filename );
      pgmwrite( filename, buf, M, N );

   }
  
   /* sum the time of all the processors to store it in total */ 
   MPI_Reduce( &endtim, &total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

   /* print the total time in rank 0*/
   if ( rank == 0 )
   {
      total = total / size;
      printf( "Elapse time for parallelising was: %.20lf\n", total );
   }

   /* shut down MPI environment: MPI_COMM_WORLD */
   MPI_Finalize();
   return 0;

}
/* boundary swaps between neighbouring processes */ 
double boundaryval( int i, int m )
{
  double val;

  val = 2.0*((double)(i-1))/((double)(m-1));
  if (i >= m/2+1) val = 2.0-val;

  return val;
}
