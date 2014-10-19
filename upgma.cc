/**
    BEAR: BEst-Aigned Rotations
    Copyright (C) 2014 Solon P. Pissis. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <sys/time.h>
#include <float.h>
#include <limits.h>

#include "beardefs.h"

#define max(a,b) ((a) > (b)) ? (a) : (b)
#define min(a,b) ((a) < (b)) ? (a) : (b)

unsigned int upgma_sim ( TPOcc ** D, unsigned int d, struct TSwitch sw, int * R )
{
	/* This is the original distanc matrix */
	double ** M;
        M = ( double ** ) calloc ( d , sizeof ( double * ) );
        if( ( M == NULL) )
        {
        	fprintf( stderr, " Error: Cannot allocate memory!\n" );
                exit( EXIT_FAILURE );
        }

	/* This is the current working distanc matrix */
	double ** N;
	N = ( double ** ) calloc ( d , sizeof ( double * ) );
	if( ( N == NULL) )
	{
		fprintf( stderr, " Error: Cannot allocate memory!\n" );
		exit( EXIT_FAILURE );
	}

	/* This is the index from N to next N */
	int * iN;
        iN = ( int * ) calloc ( d , sizeof ( int ) );
        if( ( iN == NULL) )
        {
        	fprintf( stderr, " Error: Cannot allocate memory!\n" );
                exit( EXIT_FAILURE );
        }

	/* These are the clusters: an index from M to N */
	int ** C = NULL;
        C = ( int ** ) calloc ( d , sizeof ( int * ) );
        if( ( C == NULL) )
        {
        	fprintf( stderr, " Error: Cannot allocate memory!\n" );
                exit( EXIT_FAILURE );
        }

	/* This is the number of elements per Cluster */
	int * nE;
        nE = ( int * ) calloc ( d , sizeof ( int ) );
        if( ( nE== NULL) )
        {
        	fprintf( stderr, " Error: Cannot allocate memory!\n" );
                exit( EXIT_FAILURE );
        }

	int ** pC = ( int ** ) malloc ( d * sizeof ( int * ) );	
	for ( int i = 0; i < d; i ++ )
	{
		pC[i] = NULL;
	}

	/* Initialisation */
	for ( int i = 0; i < d; i ++ )
	{
		M[i] = ( double * ) calloc ( d , sizeof ( double ) );
		if( ( M[i] == NULL) )
		{
			fprintf( stderr, " Error: Cannot allocate memory!\n" );
			exit( EXIT_FAILURE );
		}

		N[i] = ( double * ) calloc ( d , sizeof ( double ) );
		if( ( N[i] == NULL) )
		{
			fprintf( stderr, " Error: Cannot allocate memory!\n" );
			exit( EXIT_FAILURE );
		}

		C[i] = NULL;
		nE[i] = 1;
		C[i] = ( int * ) realloc ( C[i],   ( nE[i] ) * sizeof ( int ) );
		if( ( C[i] == NULL) )
		{
			fprintf( stderr, " Error: Cannot allocate memory!\n" );
			exit( EXIT_FAILURE );
		}
		C[i][0] = i;
		iN[i] = i;

		/* Initialise the arrays */
		for ( int j = 0; j < i; j++ )
		{	
			M[i][j] = ( double ) D[i][j] . err;
			N[i][j] = M[i][j];
		}
	}

	int dd = d;
	while ( dd > 1 )
	{
		/* Here we find the maximum value of N and store it at [imax, jmax] */
		unsigned int imax, jmax;
		double max = -DBL_MAX;
		for ( int i = 0; i < dd; i ++ )
		{
			for ( int j = 0; j < i; j ++ )
			{
				if ( N[i][j] > max )
				{
					max = N[i][j];
					imax = i;
					jmax = j;
				}
			}
		}

		unsigned int r = 0;
		for ( int i = 0; i < nE[imax]; i ++ )
		{
			unsigned int ii = C[imax][i];
			R[ii] = -1; //we initialise it to -1
			double br = sw . min_sim;
			for ( int j = 0; j < nE[jmax]; j ++ )
			{
				unsigned int jj = C[jmax][j];
				if ( D[ii][jj] . err >= br )
				{
					br = D[ii][jj] . err;
					r = D[ii][jj] . rot + R[jj];
					R[ii] = r;
				}
					
			}
			C[jmax] = ( int * ) realloc ( C[jmax],   ( nE[jmax] + 1 ) * sizeof ( int ) );
			if( ( C[jmax] == NULL ) )
			{
				fprintf( stderr, " Error: Cannot allocate memory!\n" );
				exit( EXIT_FAILURE );
			}
			C[jmax][nE[jmax]] = C[imax][i];
			nE[jmax]++;
		}

		#if 0
		for ( int j = 0; j < nE[jmax]; j ++ )
		{
			fprintf( stderr, " (%d, %d) ", C[jmax][j], R[C[jmax][j]] );
		}
		fprintf( stderr, "\n" );
		#endif

		#if 0
		/* Here we update the Clusters: we add the larger to the smaller */
		C[jmax] = ( int * ) realloc ( C[jmax],   ( nE[jmax] + nE[imax] ) * sizeof ( int ) );
		if( ( C[jmax] == NULL ) )
		{
			fprintf( stderr, " Error: Cannot allocate memory!\n" );
			exit( EXIT_FAILURE );
		}

		for ( int i = 0; i < nE[imax]; i ++ )
		{
			C[jmax][nE[jmax]] = C[imax][i];
			nE[jmax]++;
		}
		#endif

		for ( int i = 0; i < dd; i ++ )
		{
			pC[i] = C[i];
		}

		for ( int i = 0; i < dd; i ++ )
		{
			if ( i > imax )
			{
				C[i - 1] = pC[i];
				nE[i - 1] = nE[i];
				iN[i - 1] = iN[i];
			}
			else if ( i == imax )
			{
				free ( pC[i] );
				pC[i] == NULL;
				nE[i] = 0;
			}
			else
			{
				C[i] = pC[i];
			}
		}

		/* Compute the new array N */
		for ( int i = 0; i < dd - 1; i ++ )
		{
			for ( int j = 0; j < i; j ++ )
			{
				/* If the new cluster is involved we have to re-calculate the distances */
				if ( i == jmax || j == jmax )
				{
					double score = 0.0;
					for ( int ii = 0; ii < nE[i]; ii ++ )
					{
						for ( int jj = 0; jj < nE[j]; jj ++ )
						{
							unsigned int a = max ( C[i][ii], C[j][jj] );
							unsigned int b = min ( C[i][ii], C[j][jj] );
							score += M[a][b];
						}
					}
					N[i][j] = score / ( nE[j] * nE[i] );
				}
				else
				{
					unsigned int a = max ( iN[i], iN[j] );
					unsigned int b = min ( iN[i], iN[j] );
					N[i][j] = N[a][b];
				}
			}
		}
		
		--dd;

		/* Here we update the indices of N before the next join occurs */
		for ( int i = 0; i < dd; i ++ )
			if ( i >= imax )
				iN[i] -= 1;
	}

	free ( iN );
	free ( nE );
	for ( int i = 0; i < d; i ++ )
	{
		free ( M[i] );
		free ( N[i] );
	}
	free ( M );
	free ( C[0] );
	free ( C );
	free ( N );
	free ( pC );
	return ( 1 );
}

unsigned int upgma_dist ( TPOcc ** D, unsigned int d, struct TSwitch sw, int * R )
{
	/* This is the original distanc matrix */
	double ** M;
        M = ( double ** ) calloc ( d , sizeof ( double * ) );
        if( ( M == NULL) )
        {
        	fprintf( stderr, " Error: Cannot allocate memory!\n" );
                exit( EXIT_FAILURE );
        }

	/* This is the current working distanc matrix */
	double ** N;
	N = ( double ** ) calloc ( d , sizeof ( double * ) );
	if( ( N == NULL) )
	{
		fprintf( stderr, " Error: Cannot allocate memory!\n" );
		exit( EXIT_FAILURE );
	}

	/* This is the index from N to next N */
	int * iN;
        iN = ( int * ) calloc ( d , sizeof ( int ) );
        if( ( iN == NULL) )
        {
        	fprintf( stderr, " Error: Cannot allocate memory!\n" );
                exit( EXIT_FAILURE );
        }

	/* These are the clusters: an index from M to N */
	int ** C = NULL;
        C = ( int ** ) calloc ( d , sizeof ( int * ) );
        if( ( C == NULL) )
        {
        	fprintf( stderr, " Error: Cannot allocate memory!\n" );
                exit( EXIT_FAILURE );
        }

	/* This is the number of elements per Cluster */
	int * nE;
        nE = ( int * ) calloc ( d , sizeof ( int ) );
        if( ( nE== NULL) )
        {
        	fprintf( stderr, " Error: Cannot allocate memory!\n" );
                exit( EXIT_FAILURE );
        }

	int ** pC = ( int ** ) malloc ( d * sizeof ( int * ) );	
	for ( int i = 0; i < d; i ++ )
	{
		pC[i] = NULL;
	}

	/* Initialisation */
	for ( int i = 0; i < d; i ++ )
	{
		M[i] = ( double * ) calloc ( d , sizeof ( double ) );
		if( ( M[i] == NULL) )
		{
			fprintf( stderr, " Error: Cannot allocate memory!\n" );
			exit( EXIT_FAILURE );
		}

		N[i] = ( double * ) calloc ( d , sizeof ( double ) );
		if( ( N[i] == NULL) )
		{
			fprintf( stderr, " Error: Cannot allocate memory!\n" );
			exit( EXIT_FAILURE );
		}

		C[i] = NULL;
		nE[i] = 1;
		C[i] = ( int * ) realloc ( C[i],   ( nE[i] ) * sizeof ( int ) );
		if( ( C[i] == NULL) )
		{
			fprintf( stderr, " Error: Cannot allocate memory!\n" );
			exit( EXIT_FAILURE );
		}
		C[i][0] = i;
		iN[i] = i;

		/* Initialise the arrays */
		for ( int j = 0; j < i; j++ )
		{	
			M[i][j] = ( double ) D[i][j] . err;
			N[i][j] = M[i][j];
		}
	}

	int dd = d;
	while ( dd > 1 )
	{
		/* Here we find the minimum value of N and store it at [imin, jmin] */
		unsigned int imin, jmin;
		double min = DBL_MAX;
		for ( int i = 0; i < dd; i ++ )
		{
			for ( int j = 0; j < i; j ++ )
			{
				if ( N[i][j] < min )
				{
					min = N[i][j];
					imin = i;
					jmin = j;
				}
			}
		}


		unsigned int r = 0;
		for ( int i = 0; i < nE[imin]; i ++ )
		{
			unsigned int ii = C[imin][i];
			//we initialise it to -1
			R[ii] = -1; 
			unsigned int br = sw . k + 1;
			for ( int j = 0; j < nE[jmin]; j ++ )
			{
				unsigned int jj = C[jmin][j];
				if ( D[ii][jj] . err < br )
				{
					r = D[ii][jj] . rot + R[jj];
					R[ii] = r;
					br = D[ii][jj] . err;
				}
			}
			C[jmin] = ( int * ) realloc ( C[jmin],   ( nE[jmin] + 1 ) * sizeof ( int ) );
			if( ( C[jmin] == NULL ) )
			{
				fprintf( stderr, " Error: Cannot allocate memory!\n" );
				exit( EXIT_FAILURE );
			}	
			C[jmin][nE[jmin]] = C[imin][i];
			nE[jmin]++;
		}

		#if 0
		/* Here we update the Clusters: we add the larger to the smaller */
		C[jmin] = ( int * ) realloc ( C[jmin],   ( nE[jmin] + nE[imin] ) * sizeof ( int ) );
		if( ( C[jmin] == NULL ) )
		{
			fprintf( stderr, " Error: Cannot allocate memory!\n" );
			exit( EXIT_FAILURE );
		}

		for ( int i = 0; i < nE[imin]; i ++ )
		{
			C[jmin][nE[jmin]] = C[imin][i];
			nE[jmin]++;
		}
		#endif

		for ( int i = 0; i < dd; i ++ )
		{
			pC[i] = C[i];
		}

		for ( int i = 0; i < dd; i ++ )
		{
			if ( i > imin )
			{
				C[i - 1] = pC[i];
				nE[i - 1] = nE[i];
				iN[i - 1] = iN[i];
			}
			else if ( i == imin )
			{
				free ( pC[i] );
				pC[i] == NULL;
				nE[i] = 0;
			}
			else
			{
				C[i] = pC[i];
			}
		}

		/* Compute the new array N */
		for ( int i = 0; i < dd - 1; i ++ )
		{
			for ( int j = 0; j < i; j ++ )
			{
				/* If the new cluster is involved we have to re-calculate the distances */
				if ( i == jmin || j == jmin )
				{
					double score = 0.0;
					for ( int ii = 0; ii < nE[i]; ii ++ )
					{
						for ( int jj = 0; jj < nE[j]; jj ++ )
						{
							unsigned int a = max ( C[i][ii], C[j][jj] );
							unsigned int b = min ( C[i][ii], C[j][jj] );
							score += M[a][b];
						}
					}
					N[i][j] = score / ( nE[j] * nE[i] );
				}
				else
				{
					unsigned int a = max ( iN[i], iN[j] );
					unsigned int b = min ( iN[i], iN[j] );
					N[i][j] = N[a][b];
				}
			}
		}
		
		--dd;

		/* Here we update the indices of N before the next join occurs */
		for ( int i = 0; i < dd; i ++ )
			if ( i >= imin )
				iN[i] -= 1;
	}

	free ( iN );
	free ( nE );
	for ( int i = 0; i < d; i ++ )
	{
		free ( M[i] );
		free ( N[i] );
	}
	free ( M );
	free ( C[0] );
	free ( C );
	free ( N );
	free ( pC );
	return ( 1 );
}

