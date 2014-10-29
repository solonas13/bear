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
#include <stdint.h>
#include <float.h>
#include "beardefs.h"

#define max(a,b) ((a) > (b)) ? (a) : (b)
#define min(a,b) ((a) < (b)) ? (a) : (b)

/*
The Smith-Waterman algorithm with affine penalty scores
*/
unsigned int bcf_sw ( unsigned char * p, unsigned int m, unsigned  char * t, unsigned int n, struct TSwitch sw, TPOcc * M )
{

 	int i;
    	int j;
    	double ** D;
    	double ** I;
    	double ** T;
    	double g = sw . O;
    	double h = sw . E;
    	double matching_score = 0;

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

    	for( i = 1; i < m + 1; i++ )
    	{
		double max_score = 0;
        	for( j = 1; j < n + 1; j++ )
        	{
			D[i][j] = max ( D[i - 1][j] + h, T[i - 1][j] + g );
            		double u = D[i][j];
            		I[i][j] = max ( I[i][j - 1] + h, T[i][j - 1] + g );
            		double v = I[i][j];
            		matching_score = nuc_delta( t[j - 1], p[i - 1] );
            		if ( matching_score == ERR )
				return 0;
	            	double w = T[i - 1][j - 1] + matching_score;

			T[i][j] = max ( 0, max ( w, max ( u, v ) ) );

			if ( T[i][j] > max_score )
			{
				max_score = T[i][j];
				M[i - 1] . err  = max_score;
				M[i - 1] . rot = j - 1;
			}
        	}

    	}

        for ( i = 0; i < m + 1; i ++ )
	{
		free ( D[i] );
		free ( I[i] );
		free ( T[i] );
	}
	free ( I );
	free ( D );
	free ( T );

        return  ( 1 );
}

/*
The Smith-Waterman algorithm with affine penalty scores in linear space
*/
unsigned int bcf_sw_ls ( unsigned char * p, unsigned int m, unsigned  char * t, unsigned int n, struct TSwitch sw, TPOcc * M )
{

 	int i;
    	int j;
        double * d0;
        double * d1;
        double * t0;
        double * t1;
        double * in;
    	double g = sw . O;
    	double h = sw . E;
    	double matching_score = 0;
        
        if ( ( d0 = ( double * ) calloc ( ( n + 1 ) , sizeof ( double ) ) ) == NULL )
        {
            fprintf( stderr, " Error: 'd0' could not be allocated!\n");
            return ( 0 );
        }
        if ( ( d1 = ( double * ) calloc ( ( n + 1 ) , sizeof ( double ) ) ) == NULL  )
        {
            fprintf( stderr, " Error: 'd1' could not be allocated!\n");
            return ( 0 );
        }
        if ( ( t0 = ( double * ) calloc ( ( n + 1 ) , sizeof ( double ) ) ) == NULL )
        {
            fprintf( stderr, " Error: 't0' could not be allocated!\n");
            return ( 0 );
        }
        if ( ( t1 = ( double * ) calloc ( ( n + 1 ) , sizeof ( double ) ) ) == NULL )
        {
            fprintf( stderr, " Error: 't1' could not be allocated!\n");
            return ( 0 );
        }
        if ( ( in = ( double * ) calloc ( ( n + 1 ) , sizeof ( double ) ) ) == NULL )
        {
            fprintf( stderr, " Error: 'in' could not be allocated!\n");
            return ( 0 );
        }
        
    	for( i = 1; i < m + 1; i++ )
    	{
                double max_score = 0;
        	for( j = 1; j < n + 1; j++ )
        	{
                    double u, v, w;
                    
                    switch ( i % 2 ) {
                        
                        case 0:
                    
                            d0[j] = max ( d1[j] + h, t1[j] + g );
                            u = d0[j];
                            
                            in[j] = max ( in[j - 1] + h, t0[j - 1] + g ); //i0
                            v = in[j];
                            
                            matching_score = nuc_delta( t[j - 1], p[i - 1] );
                            if ( matching_score == ERR )
                                    return 0;
                            w = t1[j - 1] + matching_score;

                            t0[j] = max ( 0, max ( w, max ( u, v ) ) );
                            
                            if ( t0[j] > max_score )
                            {
                                    max_score = t0[j];
                                    M[i - 1] . err  = max_score;
                                    M[i - 1] . rot = j - 1;
                                    
                            }
                            
                            break;
                        
                        case 1:
                            
                            d1[j] = max ( d0[j] + h, t0[j] + g );
                            u = d1[j];
                            
                            in[j] = max ( in[j - 1] + h, t1[j - 1] + g ); //i1
                            v = in[j];
                            
                            matching_score = nuc_delta( t[j - 1], p[i - 1] );
                            if ( matching_score == ERR )
                                    return 0;
                            w = t0[j - 1] + matching_score;

                            t1[j] = max ( 0, max ( w, max ( u, v ) ) );
                            
                            if ( t1[j] > max_score )
                            {
                                    max_score = t1[j];
                                    M[i - 1] . err  = max_score;
                                    M[i - 1] . rot = j - 1;
                                    
                            }                   
                            
                            break;
                            
                    }
        	}
    	}
        
        free( d0 );
        free( d1 );
        free( t0 );
        free( t1 );
        free( in );

        return  ( 1 );
}
