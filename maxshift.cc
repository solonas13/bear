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
#include "beardefs.h"

#define max(a,b) ((a) > (b)) ? (a) : (b)
#define min(a,b) ((a) < (b)) ? (a) : (b)

/*
returns the number of 1's in an integer x
*/
inline unsigned int popcount ( WORD x )
{
	//__builtin_popcount: 32-bit
	//__builtin_popcountl: 128-bit
	//http://www.dalkescientific.com/writings/diary/archive/2011/11/02/faster_popcount_update.html
	return __builtin_popcountl( x );
}

/*
given integers a, b, c this function returns one of the integers a, b, c
with the property that it has the least number of 1's (bits set on). If there is
a draw then it returns the maximum of the two when viewed as decimal integers
*/
inline WORD bitminmax ( WORD a, WORD b, WORD c )
{
	unsigned int x , y , z; 
	WORD minimum, maximum;
	x = popcount ( a );
	y = popcount ( b );
	z = popcount ( c );

	minimum = min( x , y );
	if( z < minimum )
	{
	return c;
	}
	else if ( z == minimum && ( x != y ) )
	{
	if( x < y )
	maximum = max ( c , a );
	else
	maximum = max ( c , b );
	return maximum;
	}
	else if ( ( z == x ) && ( x == y ) )
	{
	maximum = max ( b , c );
	maximum = max ( a , maximum );
	return maximum;
	}
	else if ( z > minimum && ( x != y ) )
	{
	if( x < y )
	return a;
	else
	return b;
	}
	else /*if ( z > minimum && ( x == y ) )*/
	{
	maximum = max ( a , b );
	return maximum;
	}
}
/*
moves the bits one position to the left and enters zeros from the right
*/
inline WORD shift ( WORD a )
{
	return ( ( WORD ) a << 1 );
}
/*
shifts and truncates the leftmost bit
*/
inline WORD shiftc( WORD a, WORD x )
{
	return shift ( a & x );
}
/*
returns the hamming distance for two characters a and b
*/
inline unsigned int delta( char a, char b )
{
	if	( a == b )	return 0;
	else	return 1;
}

/*
the dynamic programming algorithm under the hamming distance model
*/
unsigned int bcf_maxshift_hd ( unsigned char * p, unsigned int m, unsigned  char * t, unsigned int n, struct TSwitch sw, TPOcc * D )
{
	WORD y; 
	WORD ** M; 		
	unsigned int i;
	unsigned int j;
	unsigned int min_err;

	if ( ( M = ( WORD ** ) calloc ( ( m + 1 ) , sizeof( WORD * ) ) ) == NULL )
	{
		fprintf( stderr, " Error: M could not be allocated!\n");
		return ( 0 );
	}
	for ( i = 0; i < m + 1; i ++ ) 
	{
		if ( ( M[i] = ( WORD * ) calloc ( ( n + 1 ) , sizeof( WORD ) ) ) == NULL )
		{
			fprintf( stderr, " Error: M could not be allocated!\n");
			return ( 0 );
		}
	}  

	y = ( ( WORD ) ( 1 ) << ( sw . w - 1 ) ) - 1;

	for ( i = 0; i < m + 1; i ++ ) 
	{  
		min_err = sw . k + 1;

		for( j = 0; j < n + 1 ; j++ )			
		{
			if( i == 0 ) continue;

			if ( j == 0 )
				M[i][j] = ( ( WORD ) 2 << ( min ( i , sw . w ) - 1 ) ) - 1;
			else
				M[i][j] = shiftc ( M[i - 1][j - 1], y ) | delta ( p[i - 1], t[j - 1] );

			if ( i >= sw . w  && j >= sw . w )
			{
				unsigned int err = popcount ( M[i][j] );

				if ( err < min_err )	
				{
					min_err = err;
					D[i - 1]  . err  = min_err;
					D[i - 1]  . rot  = j - 1;
				}
			}

                }
        }

	for ( i = 0; i < m + 1; i ++ ) 
	{
		free ( M[i] );
	}
	free ( M );	  

        return  ( 1 );
}

/*
the dynamic programming algorithm under the hamming distance model
*/
unsigned int bcf_maxshift_hd_ls ( unsigned char * p, unsigned int m, unsigned  char * t, unsigned int n, unsigned int l, unsigned int * ii, unsigned int * jj, unsigned int * distance )
{
	WORD y; 
	WORD * M0; 		
	WORD * M1; 		
	unsigned int i;
	unsigned int j;

	if ( ( M0 = ( WORD * ) calloc ( ( n + 1 ) , sizeof( WORD ) ) ) == NULL )
	{
		fprintf( stderr, " Error: M0 could not be allocated!\n");
		return ( 0 );
	}
	if ( ( M1 = ( WORD * ) calloc ( ( n + 1 ) , sizeof( WORD ) ) ) == NULL )
	{
		fprintf( stderr, " Error: M1 could not be allocated!\n");
		return ( 0 );
	}

	y = ( ( WORD ) ( 1 ) << ( l - 1 ) ) - 1;

	for ( i = 0; i < m + 1; i ++ ) 
	{  
		for( j = 0; j < n + 1 ; j++ )			
		{
			if( i == 0 ) continue;

			WORD val;

			switch ( i % 2 )
			{
				case 0 :

				if ( j == 0 )
					M0[j] = ( ( WORD ) 2 << ( min ( i , l ) - 1 ) ) - 1;
				else
					M0[j] = shiftc ( M1[j - 1], y ) | delta ( p[i - 1], t[j - 1] );

				val = M0[j];

				break;

				case 1 :

				if ( j == 0 )
					M1[j] = ( ( WORD ) 2 << ( min ( i , l ) - 1 ) ) - 1;
				else
					M1[j] = shiftc ( M0[j - 1], y ) | delta ( p[i - 1], t[j - 1] );

				val = M1[j];

				break;

				default:

				fprintf ( stderr, " Error: this should never happen!\n");

				break;	
			} 

			if ( i >= l  && j >= l )
			{
				unsigned int err = popcount ( val );

				if ( err < ( * distance ) )	
				{
					( * distance ) = err;
					( * ii )  = i - 1;
					( * jj )  = j - 1;
				}
			}

                }
        }

	free ( M0 );
	free ( M1 );	  

        return  ( 1 );
}

/*
the dynamic programming algorithm under the edit distance model
*/
unsigned int bcf_maxshift_ed_ls ( unsigned char * p, unsigned int m, unsigned  char * t, unsigned int n, unsigned int l, unsigned int * ii, unsigned int * jj, unsigned int * distance )
{
	WORD y; 
	WORD * M0; 		
	WORD * M1; 		
	unsigned int i;
	unsigned int j;

	if ( ( M0 = ( WORD * ) calloc ( ( n + 1 ) , sizeof( WORD ) ) ) == NULL )
	{
		fprintf( stderr, " Error: M0 could not be allocated!\n");
		return ( 0 );
	}
	if ( ( M1 = ( WORD * ) calloc ( ( n + 1 ) , sizeof( WORD ) ) ) == NULL )
	{
		fprintf( stderr, " Error: M1 could not be allocated!\n");
		return ( 0 );
	}
	
	y = ( ( WORD ) ( 1 ) << ( l - 1 ) ) - 1;

	for ( i = 0; i < m + 1; i ++ ) 
	{  
		for( j = 0; j < n + 1 ; j++ )			
		{
			if( i == 0 ) continue;

			WORD val;

			switch ( i % 2 )
			{
				case 0 :

				if ( j == 0 )
					M0[j] = ( ( WORD ) ( 2 ) << ( min ( i , l ) - 1 ) ) - 1;
				else if ( i <= l )
					M0[j] = bitminmax ( shift ( M1[j] ) | 1, shift( M0[j - 1] ) | 1, shift( M1[j - 1] ) | delta ( p[i - 1], t[j - 1] ) );
				else
					M0[j] = bitminmax ( shiftc( M1[j] , y ) | 1, shift ( M0[j - 1] ) | 1, shiftc ( M1[j - 1], y ) | delta ( p[i - 1], t[j - 1] ) );

				val = M0[j];

				break;

				case 1 :

				if ( j == 0 )
					M1[j] = ( ( WORD ) ( 2 ) << ( min ( i , l ) - 1 ) ) - 1;
				else if ( i <= l )
					M1[j] = bitminmax ( shift ( M0[j] ) | 1, shift( M1[j - 1] ) | 1, shift( M0[j - 1] ) | delta ( p[i - 1], t[j - 1] ) );
				else
					M1[j] = bitminmax ( shiftc( M0[j] , y ) | 1, shift ( M1[j - 1] ) | 1, shiftc ( M0[j - 1], y ) | delta ( p[i - 1], t[j - 1] ) );

				val = M1[j];

				break;

				default:

				fprintf ( stderr, " Error: this should never happen!\n");

				break;	
			} 

			if ( i >= l )
			{
				unsigned int err = popcount ( val );
				if ( err < ( * distance ) )	
				{
					( * distance ) = err;
					( * ii )  = i - 1;
					( * jj )  = j - 1;
				}
			}

                }
        }

	free ( M0 );
	free ( M1 );	  

        return  ( 1 );
}

/*
the dynamic programming algorithm under the edit distance model
*/
unsigned int bcf_maxshift_ed ( unsigned char * p, unsigned int m, unsigned  char * t, unsigned int n, struct TSwitch sw, TPOcc * D )
{
	WORD y; 
	WORD ** M; 		
	unsigned int i;
	unsigned int j;
	unsigned int min_err;

	if ( ( M = ( WORD ** ) calloc ( ( m + 1 ) , sizeof( WORD * ) ) ) == NULL )
	{
		fprintf( stderr, " Error: M could not be allocated!\n");
		return ( 0 );
	}
	for ( i = 0; i < m + 1; i ++ ) 
	{
		if ( ( M[i] = ( WORD * ) calloc ( ( n + 1 ) , sizeof( WORD ) ) ) == NULL )
		{
			fprintf( stderr, " Error: M could not be allocated!\n");
			return ( 0 );
		}
	}  

	y = ( ( WORD ) ( 1 ) << ( sw . w - 1 ) ) - 1;

	for ( i = 0; i < m + 1; i ++ ) 
	{  
		min_err = sw . k + 1;

		for( j = 0; j < n + 1 ; j++ )			
		{
			if( i == 0 ) continue;

			if ( j == 0 )
				M[i][j] = ( ( WORD ) 2 << ( min ( i , sw . w ) - 1 ) ) - 1;
			else if ( i < sw . w )
				M[i][j] = bitminmax ( shift ( M[i - 1][j] ) | 1, shift( M[i][j - 1] ) | 1, shift( M[i - 1][j - 1] ) | delta ( p[i - 1], t[j - 1] ) );
			else
				M[i][j] = bitminmax ( shiftc( M[i - 1][j] , y ) | 1, shift ( M[i][j - 1] ) | 1, shiftc ( M[i - 1][j - 1], y ) | delta ( p[i - 1], t[j - 1] ) );

			if ( i >= sw . w  )
			{
				unsigned int err = popcount ( M[i][j] );

				if ( err < min_err )	
				{
					min_err = err;
					D[i - 1]  . err  = min_err;
					D[i - 1]  . rot  = j - 1;
				}
			}

                }
        }

	for ( i = 0; i < m + 1; i ++ ) 
	{
		free ( M[i] );
	}
	free ( M );	  

        return  ( 1 );
}
