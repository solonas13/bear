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

unsigned int nw_algorithm ( unsigned char * p, unsigned int m, unsigned char * t, unsigned int n, struct TSwitch sw, double ** T )
{
	int i;
	int j;
        double **       D;
        double **       I;
	double g = sw . O;
        double h = sw . E;
	double matching_score = 0;

	D = ( double ** ) calloc ( m + 1,  sizeof ( double * ) );
	for ( int j = 0; j < m + 1; j ++ )
		D[j] = ( double * ) calloc ( n + 1,  sizeof ( double ) );
	I = ( double ** ) calloc ( m + 1,  sizeof ( int * ) );
	for ( int j = 0; j < m + 1; j ++ )
		I[j] = ( double * ) calloc ( n + 1,  sizeof ( double ) );

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

	if ( m > 0 )
		T[1][0] = g; 

        for ( i = 2; i < m + 1; i++ )
		T[i][0] = T[i - 1][0] + h;

	if ( n > 0 )
		T[0][1] = g;

        for ( j = 2; j < n + 1; j++ )
		T[0][j] = T[0][j - 1] + h;

        for( j = 1; j < n + 1; j ++ )
        {
                for( i = 1; i < m + 1; i++ )
                {

			D[i][j] = max ( D[i - 1][j] + h, T[i - 1][j] + g );
			double u = D[i][j];
			I[i][j] = max ( I[i][j - 1] + h, T[i][j - 1] + g );
			double v = I[i][j];
			matching_score = ( sw . matrix ? pro_delta( t[j - 1], p[i - 1] ) : nuc_delta( t[j - 1], p[i - 1] ) ) ;
                        if ( matching_score == ERR )
                                return 0;
			double w = T[i - 1][j - 1] + matching_score;

			T[i][j] = max ( w, max ( u, v ) );
		}
	}

	for ( j = 0; j < m + 1; j ++ )
	{
		free ( D[j] );
		free ( I[j] );
	}
	free ( D );
	free ( I );

	return ( 1 );
}


unsigned int nw_algorithm_wbt ( unsigned char * p, unsigned int m, unsigned char * t, unsigned int n, struct TSwitch sw, double ** T, int ** H )
{

	int i;
	int j;
        double **       D;
        double **       I;
	double g = sw . O;
        double h = sw . E;
	double matching_score = 0;

	D = ( double ** ) calloc ( m + 1,  sizeof ( double * ) );
	for ( int j = 0; j < m + 1; j ++ )
		D[j] = ( double * ) calloc ( n + 1,  sizeof ( double ) );
	I = ( double ** ) calloc ( m + 1,  sizeof ( int * ) );
	for ( int j = 0; j < m + 1; j ++ )
		I[j] = ( double * ) calloc ( n + 1,  sizeof ( double ) );

        for ( i = 0; i < m + 1; i++ )
	{	
		D[i][0] = -DBL_MAX;
		I[i][0] = -DBL_MAX;
	}

        for ( j = 0; j < n + 1; j++ )
	{	
		D[0][j] = -DBL_MAX;
		I[0][j] = -DBL_MAX;
	}

	for ( i = 1; i < m + 1; i++ )
    	{
      		H[i][0] = -1;
    	}
	for ( j = 1; j < n + 1; j++ )
    	{
      		H[0][j] = 1;
    	}

	T[0][0] = 0; 

	if ( m > 0 )
		T[1][0] = g; 

        for ( i = 2; i < m + 1; i++ )
		T[i][0] = T[i - 1][0] + h;

	if ( n > 0 )
		T[0][1] = g;

        for ( j = 2; j < n + 1; j++ )
		T[0][j] = T[0][j - 1] + h;

        for( j = 1; j < n + 1; j++)
        {
                for( i = 1; i < m + 1; i++ )
                {
			D[i][j] = max ( ( double ) D[i - 1][j] + h, ( double ) T[i - 1][j] + g );
			double u = D[i][j];
			I[i][j] = max ( ( double ) I[i][j - 1] + h, ( double ) T[i][j - 1] + g );
			double v = I[i][j];
			matching_score = ( sw . matrix ? pro_delta( t[j - 1], p[i - 1] ) : nuc_delta( t[j - 1], p[i - 1] ) ) ;
                        if ( matching_score == ERR )
                                return 0;
			double w = T[i - 1][j - 1] + matching_score;

			if( u > w )
			{
				if( v > u )
				{
					T[i][j] = v;
            				H[i][j] = -1; 		//an insertion
				}
				else
				{
					T[i][j] = u;
            				H[i][j] = 1;   		//a deletion
				}
			}
			else
			{
				if( v > w )
				{
					T[i][j] = v;
            				H[i][j] = -1; 		//an insertion
				}
				else
				{
					T[i][j] = w;
            				H[i][j] = 0;       	//a substitution
				}
			}
		}
	}

	for ( j = 0; j < m + 1; j ++ )
	{
		free ( D[j] );
		free ( I[j] );
	}
	free ( D );
	free ( I );

	return ( 1 );
}

unsigned int cyc_nw_algorithm_wbt ( unsigned char * p, unsigned int m, unsigned char * tt, unsigned int n, struct TSwitch sw, double * score, int * rot )
{
 	int i;
    	int j;
	int r;
        double * d0;
        double * d1;
        double * t0;
        double * t1;
        double * in;
    	double g = sw . O;
    	double h = sw . E;
    	double matching_score = 0;
	double max_score = -DBL_MAX;
        
	for ( r = 0; r < n; r++ )
	{
		unsigned char * t;	
		if ( ( t = ( unsigned char * ) calloc ( ( n + 1 ) , sizeof ( unsigned char ) ) ) == NULL )
		{
		    fprintf( stderr, " Error: 't' could not be allocated!\n");
		    return ( 0 );
		}

		create_rotation ( tt, r, t );

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

		for ( j = 0; j < n + 1; j++ )
		{	
			d0[j] = -DBL_MAX;
			in[j] = -DBL_MAX;
		}

		t0[0] = 0; 

		if ( m > 0 )
			t1[0] = g; 

		if ( n > 0 )
			t0[1] = g;

		for ( j = 2; j < n + 1; j++ )
			t0[j] = t0[j - 1] + h;

		for( i = 1; i < m + 1; i++ )
		{
			for( j = 0; j < n + 1; j++ )
			{
				double u, v, w;
			    
			    	switch ( i % 2 ) 
				{
				
					case 0:

					if ( j == 0 )
					{
						d0[j] = -DBL_MAX;
						in[j] = -DBL_MAX;
						if ( i >= 2 )
							t0[0] = t1[0] + h;
					}	
					else 
					{
		    
						d0[j] = max ( d1[j] + h, t1[j] + g );
						u = d0[j];
			    
						in[j] = max ( in[j - 1] + h, t0[j - 1] + g ); //i0
						v = in[j];
			    
						matching_score = ( sw . matrix ? pro_delta( t[j - 1], p[i - 1] ) : nuc_delta( t[j - 1], p[i - 1] ) ) ;
						if ( matching_score == ERR )
							return 0;
						w = t1[j - 1] + matching_score;

						t0[j] = max ( w, max ( u, v ) );
			    
						if ( i == m && j == n && t0[n] > max_score )
						{
							max_score = t0[j];
							( * score )  = max_score;
							( * rot ) = r;
				    
						}
					}
			    
					break;
				
					case 1:

					if ( j == 0 )
					{
						d1[j] = -DBL_MAX;
						in[j] = -DBL_MAX;
						if ( i >= 2 )
							t1[0] = t0[0] + h;
					}	
					else 
					{

						d1[j] = max ( d0[j] + h, t0[j] + g );
						u = d1[j];
			    
						in[j] = max ( in[j - 1] + h, t1[j - 1] + g ); //i1
						v = in[j];
			    
						matching_score = ( sw . matrix ? pro_delta( t[j - 1], p[i - 1] ) : nuc_delta( t[j - 1], p[i - 1] ) ) ;
						if ( matching_score == ERR )
							return 0;
						w = t0[j - 1] + matching_score;

						t1[j] = max ( w, max ( u, v ) );
			    
						if ( i == m && j == n && t1[n] > max_score )
						{
							max_score    = t1[j];
							( * score )  = max_score;
							( * rot )    = r;
				    
						} 
					}                   
			    
					break;
					    
			    	}
			}
		}
		
		free ( t );
		free( d0 );
		free( d1 );
		free( t0 );
		free( t1 );
		free( in );
	}

}
