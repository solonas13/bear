/**
    CSC: Circular Sequence Comparison
    Copyright (C) 2015 Solon P. Pissis. 

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
#include <float.h>
#include <getopt.h>
#include <assert.h>
#include <sys/time.h>
#include <divsufsort.h>
#include "globals.h"
#include "EDNAFULL.h"
#include "EBLOSUM62.h"
#include "beardefs.h"

#include <sdsl/bit_vectors.hpp>					  // include header for bit vectors

using namespace sdsl;
using namespace std;

unsigned int LCParray ( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP )
{										
	INT i=0, j=0;

	LCP[0] = 0;
	for ( i = 0; i < n; i++ ) // compute LCP[ISA[i]]
		if ( ISA[i] != 0 ) 
		{
			if ( i == 0) j = 0;
			else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
			while ( text[i+j] == text[SA[ISA[i]-1]+j] )
				j++;
			LCP[ISA[i]] = j;
		}

	return ( 1 );
}

unsigned int circular_sequence_comparison (  unsigned char * x, unsigned char * y, struct TSwitch  sw, unsigned int * rotation, unsigned int * distance )
{
	INT * SA;
	INT * LCP;
	INT * invSA;
	INT m = strlen ( ( char * ) x );
	INT n = strlen ( ( char * ) y );
	INT mmn = m + m + n;

	unsigned char * xxy;
        xxy = ( unsigned char * ) calloc( ( mmn + 1 ) , sizeof( unsigned char ) );

	strcat ( ( char * ) xxy, ( char * ) x );
	xxy[m] = '\0';
	strcat ( ( char * ) xxy, ( char * ) x );
	xxy[m + m] = '\0';
	strcat ( ( char * ) xxy, ( char * ) y );
	xxy[m + m + n] = '\0';
        
	//fprintf(stderr, " %s.\n", xxy );

        /* Compute the suffix array */
        SA = ( INT * ) malloc( ( mmn ) * sizeof( INT ) );
        if( ( SA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
                return ( 0 );
        }

        if( divsufsort( xxy, SA,  mmn ) != 0 )
        {
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
        }

        /* Compute the inverse SA array */
        invSA = ( INT * ) calloc( mmn , sizeof( INT ) );
        if( ( invSA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
                return ( 0 );
        }

        for ( INT i = 0; i < mmn; i ++ )
        {
                invSA [SA[i]] = i;
        }

	LCP = ( INT * ) calloc  ( mmn, sizeof( INT ) );
        if( ( LCP == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
                return ( 0 );
        }

        /* Compute the LCP array */
        if( LCParray( xxy, mmn, SA, invSA, LCP ) != 1 )
        {
                fprintf(stderr, " Error: LCP computation failed.\n" );
                exit( EXIT_FAILURE );
        }

	/* Ranking of q-grams and creation of x' and y' */
	int b = (int) ( m / sw . b );
	int q = sw . q;

	//fprintf(stderr, " %d %d.\n", b, q ); getchar();
	int sigma = 0;
        INT mm = m + m - q + 1; 
        INT nn = n - q + 1; 
	INT * xp; // x'
	INT * yp; // y'
        xp = ( INT * ) calloc( ( mm ) , sizeof( INT ) );
        yp = ( INT * ) calloc( ( nn ) , sizeof( INT ) );
	
	/* Here we rank the first q-gram in the suffix array */	
	if ( SA[0] >= 0 && SA[0] <= m + m  - q ) // i belongs to xx
	{
		xp[SA[0]] = sigma;
	}
	if ( SA[0] >= 2*m && SA[0] <= mmn - q ) // i belongs to y
	{
		yp[SA[0] - 2 * m] = sigma;
	}
        //fprintf(stderr, " SA[0]: %d.\n", SA[0] );

	/* Loop through the LCP array to rank the rest q-grams in the suffix array */	
	for ( INT i = 1; i < mmn; i++ )
	{
		INT lcp = LCP[i];
		INT ii = SA[i];

		//TODO: This could be optimised to ensure that sigma \in [0, m + 1]
		if ( ( ii >= 0 && ii <= m + m  - q ) || ( ii >= 2 * m && ii <= mmn - q ) )
			if ( lcp < q ) 
			{
				sigma++;
				//fprintf(stderr, " SA[i]: %d.\n", SA[i] );
			}

		if ( ii >= 0 && ii <= m + m  - q ) // i belongs to xx
		{
			xp[ii] = sigma;
		}
		if ( ii >= 2 * m && ii <= mmn - q ) // i belongs to y
		{
			yp[ii - 2 * m] = sigma;
		}
	}
        //fprintf(stderr, " sigma: %d.\n", sigma );

	#if 0
	for ( INT i = 0; i < mm; i++ )
	{
		fprintf ( stderr, "%d ", xp[i] );
	}
	fprintf ( stderr, "\n" );
	for ( INT i = 0; i < nn; i++ )
	{
		fprintf ( stderr, "%d ", yp[i] );
	}
	fprintf ( stderr, "\n" );
	#endif

	/* Partitioning x' and y' as evenly as possible */
	INT * xind; 					//this is the starting position of the fragment
	INT * xmf; 					//this is the number of q-grams in the fragment
	xind = ( INT * ) calloc ( b, sizeof ( INT ) );
	xmf = ( INT * ) calloc ( b, sizeof ( INT ) );

	for ( INT j = 0; j < b; j++ )	partitioning ( 0, j, b, m - q + 1, xmf, xind );

	INT * yind; 					//this is the starting position of the fragment
	INT * ymf; 					//this is the number of q-grams in the fragment
	yind = ( INT * ) calloc ( b, sizeof ( INT ) );
	ymf = ( INT * ) calloc ( b, sizeof ( INT ) );

	for ( INT j = 0; j < b; j++ )	partitioning ( 0, j, b, nn, ymf, yind );

	#if 0
	for ( INT i = 0; i < b; i++ )
	{
		fprintf ( stderr, "(%d %d) ", xind[i], xmf[i] );
	}
	fprintf ( stderr, "\n" );
	#endif

	/* Allocate the diff vector */
	INT ** diff;
	diff = ( INT ** ) calloc ( b, sizeof ( INT * ) );
	for ( INT i = 0; i < b; i++ )	diff[i] = ( INT * ) calloc ( sigma + 1, sizeof ( INT ) );

	/* Step 1: Create diff, pvy, and D_0 */
	INT * D;
	D = ( INT * ) calloc ( b, sizeof ( INT * ) );
	for ( INT i = 0; i < b; i++ )
	{
		for ( INT j = yind[i]; j < yind[i] + ymf[i]; j++ )
		{
			diff[i][yp[j]]++;	
			D[i]++;
		}
	}	
	//fprintf ( stderr, "D0 = %d\n", D[0] ); getchar();

	/* Step 2: Compute the distances for position 0 */
	int min_dist = 0;
	for ( INT i = 0; i < b; i++ )	//first window
	{	
		for ( INT j = xind[i]; j < xind[i] + xmf[i]; j++ )
		{
			diff[i][xp[j]]--;
			if ( diff[i][xp[j]] >= 0 )
			{
				D[i]--;
			}
			else
			{
				D[i]++;
			}
		}
		min_dist += D[i];
	}
	//fprintf ( stderr, "D0 = %d\n", D[0] ); getchar();

	/* Step 3: Compute the rest of the distances */
	int rot = 0;
	for ( INT i = 1; i < m; i++ )	//all the rest windows
	{
		int dist = 0;
		for ( INT j = 0; j < b; j++ )
		{
			diff[j][xp[i - 1 + xind[j]]]++; //letter out
			diff[j][xp[i - 1 + xind[j] + xmf[j]]]--; //letter in

			//For the letter we take out
			if ( diff[j][xp[i - 1 + xind[j]]] <= 0 )	
			{
				D[j]--;
			}
			else
			{
				D[j]++;
			}

			//For the letter we add in
			if ( diff[j][xp[i - 1 + xind[j] + xmf[j]]] < 0 )	
			{
				D[j]++;
			}
			else
			{							
				D[j]--;
			}
			dist += D[j];
		}
		//fprintf ( stderr, "dist = %d\n", dist );
		if ( dist < min_dist )
		{
			rot = i;
			min_dist = dist;
		}
	}
	( * distance ) = ( unsigned int ) min_dist; 
	( * rotation ) = ( unsigned int ) rot;

	if ( sw . P > 0 )
	{
		sacsc_refinement ( x, y, sw, rotation, distance );
	}

	/* De-allocate the memory */	
	free ( D );
	free ( xp );
	free ( yp );
	free ( xind );
	free ( xmf );
	free ( yind );
	free ( ymf );
	free ( xxy );
	free ( invSA );
	free ( SA );
	free ( LCP );
	for ( INT i = 0; i < b; i++ )	free ( diff[i] );
	free ( diff );
	return ( 1 );
}

void partitioning ( INT i, INT j, INT f, INT m, INT * mf, INT * ind )
{
    	INT modulo = m % f;
    	double nf = ( double ) ( m ) / f;
    	INT first;
    	INT last;
	if ( j < modulo )
	{
		first = j * ( ceil( nf ) );
		last = first + ceil( nf ) - 1;
	}
	else
	{
		first = j * ( floor( nf ) ) + modulo;
		last = first + floor( nf ) - 1;
	}
	ind[j + i * f] = first;
	mf[j + i * f] = last - first + 1;
}

/**
 * Needleman-Wunsch for the sacsc_refinement method
 */
unsigned int nw_refine ( unsigned char * p, unsigned int m, unsigned char * t, unsigned int n, double * score, struct TSwitch sw )
{
	int i, j;
	double g = sw . O;
	double h = sw . E;
	double u, v, w;

	double ** T;
	double ** I;
	double ** D;

	if ( ( T = ( double ** ) calloc ( ( m + 1 ) , sizeof( double * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: T could not be allocated!\n");
                return ( 0 );
        }
        for ( i = 0; i < m + 1; i ++ )
        {
                if ( ( T[i] = ( double * ) calloc ( ( n + 1 ) , sizeof( double ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: T could not be allocated!\n");
                        return ( 0 );
                }
        }
	if ( ( I = ( double ** ) calloc ( ( m + 1 ) , sizeof( double * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: I could not be allocated!\n");
                return ( 0 );
        }
        for ( i = 0; i < m + 1; i ++ )
        {
                if ( ( I[i] = ( double * ) calloc ( ( n + 1 ) , sizeof( double ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: I could not be allocated!\n");
                        return ( 0 );
                }
        }
	if ( ( D = ( double ** ) calloc ( ( m + 1 ) , sizeof( double * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: D could not be allocated!\n");
                return ( 0 );
        }
        for ( i = 0; i < m + 1; i ++ )
        {
                if ( ( D[i] = ( double * ) calloc ( ( n + 1 ) , sizeof( double ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: D could not be allocated!\n");
                        return ( 0 );
                }
        }

        for ( i = 0; i < m + 1; i++ )
	{	
		D[i][0] = -DBL_MAX;
		I[i][0] = -DBL_MAX;
	}

        for ( j = 1; j < n + 1; j++ )
	{	
		D[0][j] = -DBL_MAX;
		I[0][j] = -DBL_MAX;
	}

	T[0][0] = 0;
	T[1][0] = g; 

        for ( i = 2; i < m + 1; i++ )
		T[i][0] = T[i - 1][0] + h;

	T[0][1] = g;

        for ( j = 2; j < n + 1; j++ )
		T[0][j] = T[0][j - 1] + h;

	for( i = 1; i < m + 1; i++ )
	{
        	for( j = 1; j < n + 1; j++ )
        	{
			D[i][j] = max ( D[i - 1][j] + h, T[i - 1][j] + g );
			u = D[i][j];

			I[i][j] = max ( I[i][j - 1] + h, T[i][j - 1] + g );
			v = I[i][j];

			w = T[i - 1][j - 1] + ( sw . matrix == 0 ) ? nuc_delta ( t[j - 1], p[i - 1] ) : pro_delta ( t[j - 1], p[i - 1] );

			T[i][j] = max ( w, max ( u, v ) );
        	}
    	}
	( * score ) = T[m][n];

        for ( i = 0; i < m + 1; i ++ )
	{
		free ( D[i] );
		free ( I[i] );
		free ( T[i] );
	}
	free ( I );
	free ( D );
	free ( T );
	
	return EXIT_SUCCESS;
}

/*
 * The new refine method uses 3*sw.P*sw.b long seqs and match_char (DEL)
 * represents a neutral match score regardless of the base being matched, thus
 * reducing gap opening and promoting substitution instead (see nuc_delta())
 */
unsigned int sacsc_refinement (  unsigned char * x, unsigned char * y, struct TSwitch  sw, unsigned int * rotation, unsigned int * distance )
{
	init_substitution_score_tables ();
	unsigned int rot = ( * rotation );
	unsigned int dist = ( * distance );

	unsigned int m = strlen ( ( char * ) x );
	unsigned int n = strlen ( ( char * ) y );
	unsigned char * xr;
	xr = ( unsigned char * ) calloc( ( m + 1 ) , sizeof( unsigned char ) );
	create_rotation ( x, rot, xr );
	
	unsigned char * X;
	unsigned char * Y;

	unsigned int sl = sw . P * ( sw . b ); //section length
	sl = min ( sl, min ( m/2, n/2 ) );

	X = ( unsigned char * ) calloc( ( 3 * sl + 1 ) , sizeof( unsigned char ) );
	Y = ( unsigned char * ) calloc( ( 3 * sl + 1 ) , sizeof( unsigned char ) );


	memcpy ( &X[0], &xr[0], sl );
	for ( int i = 0; i < sl; i++ )
		X[sl + i] = DEL;
	memcpy ( &X[sl + sl], &xr[m - sl], sl );
	X[3 * sl] = '\0';
	
	memcpy ( &Y[0], &y[0], sl );
	for ( int i = 0; i < sl; i++ )
		Y[sl + i] = DEL;
	memcpy ( &Y[sl + sl], &y[n - sl], sl );
	Y[3 * sl] = '\0';

	unsigned int mm = sl + sl + sl;
	unsigned int nn = sl + sl + sl;

	double score = -DBL_MAX;
	double max_score = score;
	unsigned int rrot = 0;
	unsigned char * Xr = ( unsigned char * ) calloc( ( 3 * sl + 1 ) , sizeof( unsigned char ) );
	for ( int i = 0; i < mm; i++ )
	{
		if ( i >= sl && i < 2 * sl )
			continue;
	
		Xr[0] = '\0';
		create_rotation ( X, i, Xr );

		nw_refine ( Xr, mm , Y, nn, &score, sw );
		if ( score > max_score )
		{
			max_score = score;
			rrot = i;
		}	 
	}
	free ( Xr);

	int final_rot;
	if ( rrot < sl )
	{
		final_rot = rot + rrot;
	}
	else
	{
		final_rot = rot - ( 3 * sl - rrot );
	}

	if ( final_rot > ( int ) m )
	{
		( * rotation ) = final_rot % m;	
	}
	else if ( final_rot < 0 )
	{
		( * rotation ) = m + final_rot;
	}
	else
		( * rotation ) = final_rot;

	free ( xr );
	free ( X );
	free ( Y );
	return EXIT_SUCCESS;
}
