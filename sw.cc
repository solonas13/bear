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
the Smith-Waterman algorithm
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
                        matching_score = ( sw . matrix ? pro_delta( t[j - 1], p[i - 1] ) : nuc_delta( t[j - 1], p[i - 1] ) ) ;
            		if ( matching_score == ERR )
            		{
				return 0;
            		}
	            	double w = T[i - 1][j - 1] + matching_score;

            		if( u <= 0 && v <= 0 && w <= 0 )
            		{
                		T[i][j] = 0;
            		}
            		else if ( u > w )
            		{
                		if( v > u )
                		{
                    			T[i][j] = v;
                		}
                		else
                		{
                    			T[i][j] = u;
                		}
            		}
            		else
            		{
                		if( v > w )
                		{
                    			T[i][j] = v;
                		}
				else
				{
					T[i][j] = w;
				}
            		}

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
