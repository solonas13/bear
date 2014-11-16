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
#include <float.h>
#include <limits.h>
#include <vector>
#include <tuple>
#include <algorithm>

#include "beardefs.h"
#include "filter.h"
#include "aca.h"
#include "globals.h"

#define OPT_LEN_FRA 20
#define GAPS_PEN    1
#define MIS_PEN     1

using namespace std;

/* Definition of global variables defined in globals.h */
int * gF  = NULL;
int * gP  = NULL;
int gMatches = 0;
int gMax_alloc_matches = ALLOC_SIZE;

unsigned int macsmf_ed( unsigned char ** x, unsigned char * t, struct TSwitch sw, TPOcc *** POcc, unsigned int ** NOcc )
{

	unsigned int f = 2 * sw . k + 4;	//this is the number of fragments
	unsigned int d = 0;			//this is the number of patterns in the set.
	unsigned int i = 0;
	unsigned char ** Tmp;
	for ( Tmp = x; *Tmp; Tmp++, d++ );

        unsigned int * MOcc;
        MOcc = ( unsigned int * ) calloc ( d , sizeof ( unsigned int ) );
        if( ( MOcc == NULL) )
        {
        	fprintf( stderr, " Error: Cannot allocate memory!\n" );
                exit( EXIT_FAILURE );
        }
	
	unsigned int n = strlen ( ( char * ) t );		//this is the length of the text.

	unsigned int * m; 
	m = ( unsigned int * ) calloc ( d, sizeof ( unsigned int ) );
	if ( m == NULL )
	{
		fprintf ( stderr, " Error: Cannot allocate memory for m!\n" );
		return 0;
	}

	for ( i = 0; i < d; i++ )
	{
		m[i] = strlen ( ( char * ) x[i] );		//this is the length of the pattern
	}

	unsigned char ** xx; 
	xx = ( unsigned char ** ) malloc ( d * sizeof ( unsigned char * ) );
	if ( xx == NULL )
	{
		fprintf ( stderr, " Error: Cannot allocate memory for xx!\n" );
		return 0;
	}

	for ( i = 0; i < d; i++ )
	{
		unsigned int mm = 2 * m[i] - 1;
		xx[i] = ( unsigned char * ) calloc ( mm + 1, sizeof ( unsigned char ) );
		if ( xx[i] == NULL )
		{
			fprintf ( stderr, " Error: Cannot allocate memory for xx[%d]!\n", i );
			return 0;
		}
		memmove ( &xx[i][0], x[i], m[i] );
		memmove ( &xx[i][m[i]], x[i], m[i] - 1 );
		xx[i][mm] = '\0';
		if ( ( int ) ( mm / f ) >= OPT_LEN_FRA )
			f = ( int ) ( mm / OPT_LEN_FRA );
		else
			f = 2 * sw . k + 4;	
	}

	int * ind;			//this is the starting position of the fragment
	int * mf;			//this is the length of the fragment
	int * whe;			//this is which pattern the fragment is extracted from

	ind = ( int * ) calloc ( f * d, sizeof ( int ) );
	mf = ( int * ) calloc ( f * d, sizeof ( int ) );
	whe = ( int * ) calloc ( f * d, sizeof ( int ) );
	
	for ( int i = 0; i < d; i++ )
	{
		unsigned int mm = 2 * m[i] - 1;
		for ( int j = 0; j < f; j++ )
		{
			fragments ( i, j, f, mm, mf, ind );
			whe[i * f + j] = i;
		}
	}

	/* Check whether there exist duplicated fragments */
	fprintf ( stderr, " Checking whether there exist any duplicates fragments in the patterns\n" );
	char ** seqs;
	int   * dups;		
        dups  = ( int * ) calloc ( f * d, sizeof ( int ) );
	unsigned int uniq;
	uniq = extract_dups ( xx, d, m, f, mf, ind, dups );

	int   	* d_occ;	// d_occ[i] stores the id of the next fragment that is equal to itself a la linked-list manner
	int   	* l_occ;
        l_occ = ( int * ) calloc ( f * d, sizeof ( int ) );
        d_occ = ( int * ) calloc ( f * d, sizeof ( int ) );
	for ( int i = 0; i < f * d; i++ )	
	{
		d_occ[i] = -1;
		l_occ[i] = -1;
	}

	/* In case there exist duplicated fragmnents */
	if ( uniq < f * d )
	{
               	seqs = ( char ** ) calloc ( f * d, sizeof ( char * ) );
		for ( int i = 0; i < d; i++ )	
		{
			for ( int j = 0; j < f; j++ )	
			{
				unsigned int f_id = i * f + j;

				/* Add the fragment once */
				if ( dups[f_id] < 0 )
				{
					seqs[f_id] = ( char * ) calloc ( mf[f_id] + 1, sizeof ( char ) );
					memmove ( &seqs[f_id][0], &xx[i][ ind[f_id] ], mf[f_id] );
					seqs[f_id][mf[f_id]] = '\0';
				}
				else //add nothing since it is already added 
				{
					seqs[f_id] = ( char * ) calloc ( 1, sizeof ( char ) );
					seqs[f_id][0] = '\0';

					if ( l_occ[ dups[f_id] ] == -1 )		//if it the first duplicated fragment
						d_occ[ dups[f_id] ] = f_id;
					else
						d_occ[ l_occ[ dups[f_id] ] ] = f_id;
					l_occ[ dups[f_id] ] = f_id;
				}		
			}
		}	
	}
	else //add all the fragments since there exist no duplicated fragment
	{
                seqs = ( char ** ) calloc ( f * d, sizeof ( char * ) );
		for( int i = 0; i < d; ++i )
        	{
			for ( int j = 0; j < f; j++ )	
			{
				unsigned int f_id = i * f + j;
                		seqs[f_id] = ( char * ) calloc ( mf[f_id] + 1, sizeof ( char ) );
               	 		memmove ( &seqs[f_id][0], &xx[i][ ind[f_id] ], mf[f_id] );
                		seqs[f_id][mf[f_id]] = '\0';
			}
        	}
	}

	int * F = NULL;
	F = ( int * ) realloc ( F,  ( ALLOC_SIZE ) * sizeof ( int ) );
	if ( F == NULL )
	{
		fprintf ( stderr, " Error: Cannot allocate memory for F.\n" );
		return 0;
	}

	int * P = NULL;
	P = ( int * ) realloc ( P,  ( ALLOC_SIZE ) * sizeof ( int ) );
	if ( P == NULL ) 
	{
		fprintf ( stderr, " Error: Cannot allocate memory for P.\n" );
		return 0;
	}

	int matches;

        gF  = F;
        gP  = P;

	/* Aho Corasick Automaton */
	fprintf ( stderr, " Building up the AC automaton\n" );
	filtering ( ( char * ) t, n, ( char ** ) seqs, f * d );

	F = gF;
	P = gP;
	matches = gMatches;
	
	unsigned char * tr;
	tr = ( unsigned char * ) calloc ( n + 1,  sizeof ( unsigned char ) );
	if ( tr == NULL )
	{
		fprintf ( stderr, " Error: Cannot allocate memory for Tr.\n" );
		return 0;
	}
	for ( int j = 0; j < n; j++ )
		tr[j] = t[n - j - 1];
	tr[n] = '\0';

	unsigned char ** xxr; 
	xxr = ( unsigned char ** ) malloc ( d * sizeof ( unsigned char * ) );
	if ( xxr == NULL )
	{
		fprintf ( stderr, " Error: Cannot allocate memory for xx!\n" );
		return 0;
	}

	for ( i = 0; i < d; i++ )
	{
		unsigned int mm = 2 * m[i] - 1;
		xxr[i] = ( unsigned char * ) calloc ( mm + 1, sizeof ( unsigned char ) );
		if ( xxr[i] == NULL )
		{
			fprintf ( stderr, " Error: Cannot allocate memory for xx[%d]!\n", i );
			return 0;
		}
		for ( int j = 0; j < mm; j++ )
			xxr[i][j] = xx[i][mm - j - 1];
		xxr[i][mm] = '\0';
	}


	fprintf ( stderr, " Matching starts\n" );
	for ( int i = 0; i < matches; i++ )
	{
		int jj = F[i];		// this is the ID of fragment.
		do
		{
			int w = whe[jj];	// this is which pattern the fragment is extracted from.
			unsigned int mm = m[w] * 2 - 1;

			int ii = P[i];
			int iir = n - P[i] - mf[jj];

			int r = mm - ind[jj] - mf[jj];
			int l = ind[jj];
			int exr = min( n - P[i] - mf[jj] , r + sw . k );
			int exl = min( P[i], l + sw . k );

			int pxx = ind[jj] + mf[jj];
			int pt = ii + mf[jj];
			int pxxr = mm - ind[jj];
			int ptr = iir + mf[jj];

			int * Sr;
			int * Sl;
			int * Slr;
			int * Pl;
			int * Plr;
			
			/* the right part*/
			if ( r > 0 )
			{
				Sr = ( int * ) calloc ( r,  sizeof ( int ) );

				unsigned int ** D;
				D = ( unsigned int ** ) calloc ( r + 1,  sizeof ( unsigned int * ) );
				for ( int j = 0; j < r + 1; j ++ )
					D[j] = ( unsigned int * ) calloc ( exr + 1,  sizeof ( unsigned int ) );

				edit_distance ( &xx[w][pxx], r, &t[pt], exr, sw, D );
				r_errors_vec	  ( D, r, exr, sw, Sr );

				for ( int j = 0; j < r + 1; j ++ )
				{
					free ( D[j] );
				}
				free ( D );
			}

			/* the left part */
			if ( l > 0 )
			{
				Sl = ( int * ) calloc ( l,  sizeof ( int ) );
				Pl = ( int * ) calloc ( l,  sizeof ( int ) );
				Slr = ( int * ) calloc ( l,  sizeof ( int ) );
				Plr = ( int * ) calloc ( l,  sizeof ( int ) );

				unsigned int ** D;
				int ** H;

				D = ( unsigned int ** ) calloc ( l + 1,  sizeof ( unsigned int * ) );
				for ( int j = 0; j < l + 1; j ++ )
					D[j] = ( unsigned int * ) calloc ( exl + 1,  sizeof ( unsigned int ) );

				H = ( int ** ) calloc ( l + 1,  sizeof ( int * ) );
				for ( int j = 0; j < l + 1; j ++ )
					H[j] = ( int * ) calloc ( exl + 1,  sizeof ( int ) );

				edit_distance_wbt ( &xxr[w][pxxr], l, &tr[ptr], exl, sw, D, H );
				l_errors_vec	  ( D, l, exl, H, sw, Slr, Plr );

				for ( int j = 0; j < l; j ++ )
				{
					Sl[l - j - 1] = Slr[j];
					Pl[l - j - 1] = Plr[j];
				}
				for ( int j = 0; j < l + 1; j ++ )
				{
					free ( D[j] );
					free ( H[j] );
				}
				free ( D );
				free ( H );
			}

			int * M = ( int * ) calloc ( mm, sizeof ( int ) );
			unsigned int a = 0;
			unsigned int b = 0;

			for ( int j = 0; j < mm; j++ )
			{
				if ( j < l )
					M[j] = Sl[a++];
				if ( j >= l + mf[jj] )
					M[j] = Sr[b++];
			}		

                        int left  = max ( 0, ( int ) ( ind[jj] + mf[jj] - m[w] )  );	//leftmost rotation j
                        int right = min ( ind[jj], m[w] - 1 );				//rightmost rotation j

			for ( int j = left; j <= right; j++ )
			{
				int dist = 0;
				int led = 0;
				int red = 0;
				int pl = 0;

                                if ( j < ind[jj] )
				{
                                        led = M[j];
					pl = Pl[j];
				}

                                if ( j + m[w] > ind[jj] + mf[jj] )
				{
                                        red = M[j + m[w] - 1];
				}

				dist = led + 0 + red;

				if ( dist <= sw . k )
				{
					int txt_chars = l + pl;
					int pos = ii - txt_chars + j;

					if ( ( * NOcc )[w] >= MOcc[w] )
					{
						( * POcc ) [w] = ( TPOcc * ) realloc ( ( * POcc )[w],   ( MOcc[w] + ALLOC_SIZE ) * sizeof ( TPOcc ) );
						MOcc[w] += ALLOC_SIZE;
					}

					( *POcc )[w][ ( * NOcc )[w] ] . pos = pos;	//text index
					( *POcc )[w][ ( * NOcc )[w] ] . err = dist;	
					( *POcc )[w][ ( * NOcc )[w] ] . rot = j;	
					( * NOcc )[w] = ( * NOcc )[w] + 1;
				}
			}	

			if ( r > 0 )
				free ( Sr );

			if ( l > 0 )
			{
				free ( Pl );
				free ( Sl );
				free ( Plr );
				free ( Slr );
			}

			free ( M );
		
			jj = d_occ[jj];
			
		} while ( jj != -1 );
	}

	free ( F );
	free ( P );

	/* A bit of post-processing to get rid of the duplicated occurrences */
	for ( i = 0; i < d; i++ )
	{
		vector<mytuple> data;

		int j;
		for ( j = 0; j < ( * NOcc )[i]; j++ )
  			data.push_back(make_tuple(( *POcc )[i][j] . pos, ( *POcc )[i][j] . err, ( *POcc )[i][j] . rot));

  		std::sort( data.begin(), data.end() );
		
		j = 0;
		for(vector<mytuple>::iterator iter = data.begin(); iter != data.end(); iter++)
		{
    			( *POcc )[i][j] . pos = get<0>(*iter);
			( *POcc )[i][j] . err = get<1>(*iter);
			( *POcc )[i][j] . rot = get<2>(*iter);
			j++;
  		}

		/* Here we remove the occurrences by imlpementing the C++ unique function */
		unsigned int ue = 0;
                TPOcc * Last =  unique ( &( ( *POcc )[i][0] ), &( ( *POcc )[i][( * NOcc )[i]] ) );
	
		TPOcc * Current;
		for ( Current = &( ( *POcc )[i][0] ); Current != Last; ++Current, ++ ue );
                ( * NOcc )[i] = ue;

		/* Here we resize the vector to save space */
		if ( ( ( * NOcc )[i] ) > 0  )
			( * POcc ) [i]   = ( TPOcc * ) realloc ( ( * POcc ) [i], ( ( * NOcc )[i] ) * sizeof ( TPOcc ) );
	} 

        free ( MOcc );

	free ( m );
	for ( int i = 0; i < d; i++ )
	{
		free( xx[i] );
		free( xxr[i] );
	}	
	free ( xx );
	free ( xxr );

	free ( tr );

	free ( mf );
	free ( ind );
	free ( whe );

	for ( int i = 0; i < f * d; i++ )
		free ( seqs[i] );
	free ( seqs );

	free ( l_occ );
	free ( d_occ );
	free ( dups );

	return 1;
}

unsigned int macsmf_hd( unsigned char ** x, unsigned char * t, struct TSwitch sw, TPOcc *** POcc, unsigned int ** NOcc )
{
	unsigned int f = 2 * sw . k + 4;	//this is the number of fragments
	unsigned int d = 0;			//this is the number of patterns in the set.
	unsigned int i = 0;
	unsigned char ** Tmp;
	for ( Tmp = x; *Tmp; Tmp++, d++ );

        unsigned int * MOcc;
        MOcc = ( unsigned int * ) calloc ( d , sizeof ( unsigned int ) );
        if( ( MOcc == NULL) )
        {
        	fprintf( stderr, " Error: Cannot allocate memory!\n" );
                exit( EXIT_FAILURE );
        }
	
	unsigned int n = strlen ( ( char * ) t );		//this is the length of the text.

	unsigned int * m; 
	m = ( unsigned int * ) calloc ( d, sizeof ( unsigned int ) );
	if ( m == NULL )
	{
		fprintf ( stderr, " Error: Cannot allocate memory for m!\n" );
		return 0;
	}

	for ( i = 0; i < d; i++ )
	{
		m[i] = strlen ( ( char * ) x[i] );		//this is the length of the pattern
	}

	unsigned char ** xx; 
	xx = ( unsigned char ** ) malloc ( d * sizeof ( unsigned char * ) );
	if ( xx == NULL )
	{
		fprintf ( stderr, " Error: Cannot allocate memory for xx!\n" );
		return 0;
	}

	for ( i = 0; i < d; i++ )
	{
		unsigned int mm = 2 * m[i] - 1;
		xx[i] = ( unsigned char * ) calloc ( mm + 1, sizeof ( unsigned char ) );
		if ( xx[i] == NULL )
		{
			fprintf ( stderr, " Error: Cannot allocate memory for xx[%d]!\n", i );
			return 0;
		}
		memmove ( &xx[i][0], x[i], m[i] );
		memmove ( &xx[i][m[i]], x[i], m[i] - 1 );
		xx[i][mm] = '\0';

		if ( ( int ) ( mm / f ) >= OPT_LEN_FRA )
			f = ( int ) ( mm / OPT_LEN_FRA );
		else
			f = 2 * sw . k + 4;	
	}

	int * ind;			//this is the starting position of the fragment
	int * mf;			//this is the length of the fragment
	int * whe;			//this is which pattern the fragment is extracted from

	ind = ( int * ) calloc ( f * d, sizeof ( int ) );
	mf = ( int * ) calloc ( f * d, sizeof ( int ) );
	whe = ( int * ) calloc ( f * d, sizeof ( int ) );
	
	for ( int i = 0; i < d; i++ )
	{
		unsigned int mm = 2 * m[i] - 1;
		for ( int j = 0; j < f; j++ )
		{
			fragments ( i, j, f, mm, mf, ind );
			whe[i * f + j] = i;
		}
	}

	/* Check whether there exist duplicated fragments */
	fprintf ( stderr, " Checking whether there exist any duplicates fragments in the patterns\n" );
	char ** seqs;
	int   * dups;		

        dups  = ( int * ) calloc ( f * d, sizeof ( int ) );
	unsigned int uniq;
	uniq = extract_dups ( xx, d, m, f, mf, ind, dups );

	int   	* d_occ;	// d_occ[i] stores the id of the next fragment that is equal to itself a la linked-list manner
	int   	* l_occ;
        l_occ = ( int * ) calloc ( f * d, sizeof ( int ) );
        d_occ = ( int * ) calloc ( f * d, sizeof ( int ) );
	for ( int i = 0; i < f * d; i++ )	
	{
		d_occ[i] = -1;
		l_occ[i] = -1;
	}

	/* In case there exist duplicated fragmnents */
	if ( uniq < f * d )
	{
               	seqs = ( char ** ) calloc ( f * d, sizeof ( char * ) );
		for ( int i = 0; i < d; i++ )	
		{
			for ( int j = 0; j < f; j++ )	
			{
				unsigned int f_id = i * f + j;

				/* Add the fragment once */
				if ( dups[f_id] < 0 )
				{
					seqs[f_id] = ( char * ) calloc ( mf[f_id] + 1, sizeof ( char ) );
					memmove ( &seqs[f_id][0], &xx[i][ ind[f_id] ], mf[f_id] );
					seqs[f_id][mf[f_id]] = '\0';
				}
				else //add nothing since it is already added 
				{
					seqs[f_id] = ( char * ) calloc ( 1, sizeof ( char ) );
					seqs[f_id][0] = '\0';

					if ( l_occ[ dups[f_id] ] == -1 )		//if it the first duplicated fragment
						d_occ[ dups[f_id] ] = f_id;
					else
						d_occ[ l_occ[ dups[f_id] ] ] = f_id;
					l_occ[ dups[f_id] ] = f_id;
				}		
			}
		}	
	}
	else //add all the fragments since there exist no duplicated fragment
	{
                seqs = ( char ** ) calloc ( f * d, sizeof ( char * ) );
		for( int i = 0; i < d; ++i )
        	{
			for ( int j = 0; j < f; j++ )	
			{
				unsigned int f_id = i * f + j;
                		seqs[f_id] = ( char * ) calloc ( mf[f_id] + 1, sizeof ( char ) );
               	 		memmove ( &seqs[f_id][0], &xx[i][ ind[f_id] ], mf[f_id] );
                		seqs[f_id][mf[f_id]] = '\0';
			}
        	}
	}

	int * F = NULL;
	F = ( int * ) realloc ( F,  ( ALLOC_SIZE ) * sizeof ( int ) );
	if ( F == NULL )
	{
		fprintf ( stderr, " Error: Cannot allocate memory for F.\n" );
		return 0;
	}

	int * P = NULL;
	P = ( int * ) realloc ( P,  ( ALLOC_SIZE ) * sizeof ( int ) );
	if ( P == NULL ) 
	{
		fprintf ( stderr, " Error: Cannot allocate memory for P.\n" );
		return 0;
	}

	int matches;

        gF  = F;
        gP  = P;

	/* Aho Corasick Automaton */
	fprintf ( stderr, " Building up the AC automaton\n" );
	filtering ( ( char * ) t, n, ( char ** ) seqs, f * d );

        F  = gF;
        P  = gP;
	matches = gMatches;

	unsigned char * tr;
	tr = ( unsigned char * ) calloc ( n + 1,  sizeof ( unsigned char ) );
	if ( tr == NULL )
	{
		fprintf ( stderr, " Error: Cannot allocate memory for Tr.\n" );
		return 0;
	}
	for ( int j = 0; j < n; j++ )
		tr[j] = t[n - j - 1];
	tr[n] = '\0';

	unsigned char ** xxr; 
	xxr = ( unsigned char ** ) malloc ( d * sizeof ( unsigned char * ) );
	if ( xxr == NULL )
	{
		fprintf ( stderr, " Error: Cannot allocate memory for xx!\n" );
		return 0;
	}

	for ( i = 0; i < d; i++ )
	{
		unsigned int mm = 2 * m[i] - 1;
		xxr[i] = ( unsigned char * ) calloc ( mm + 1, sizeof ( unsigned char ) );
		if ( xxr[i] == NULL )
		{
			fprintf ( stderr, " Error: Cannot allocate memory for xx[%d]!\n", i );
			return 0;
		}
		for ( int j = 0; j < mm; j++ )
			xxr[i][j] = xx[i][mm - j - 1];
		xxr[i][mm] = '\0';
	}

	fprintf ( stderr, " Matching starts\n" );
	for ( int i = 0; i < matches; i++ )
	{
		int jj = F[i];		// this is the ID of fragment.
		do
		{
			int w = whe[jj];	// this is which pattern the fragment is extracted from.
			unsigned int mm = m[w] * 2 - 1;

			int ii = P[i];
			int iir = n - P[i] - mf[jj];

			
			int r = mm - ind[jj] - mf[jj];
			int l = ind[jj];
			int exr = min( n - P[i] - mf[jj] , r );
			int exl = min( P[i], l );

			int pxx = ind[jj] + mf[jj];
			int pt = ii + mf[jj];
			int pxxr = mm - ind[jj];
			int ptr = iir + mf[jj];

			unsigned int * Sr;
			unsigned int * Sl;
			unsigned int * Slr;
			
			/* the right part*/
			if ( r > 0 )
			{
				Sr = ( unsigned int * ) calloc ( r,  sizeof ( unsigned int ) );

				/* If there is not text on the right to extend */
				if ( exr == 0 )
				{
					for ( int j = 0; j < r; j ++ )
						Sr[j] = sw . k + 1;
				}
				else
				{
					hamming_distance ( &xx[w][pxx], r, &t[pt], exr, sw , Sr );
				}
			}

			/* the left part */
			if ( l > 0 )
			{
				Sl = ( unsigned int * ) calloc ( l,  sizeof ( unsigned int ) );
				Slr = ( unsigned int * ) calloc ( l,  sizeof ( unsigned int ) );

				/* If there is not text on the left to extend */
				if ( exl == 0 )
				{
					for ( int j = 0; j < l; j ++ )
						Sl[j] = sw . k + 1;
				}
				else
				{
					hamming_distance ( &xxr[w][pxxr], l, &tr[ptr], exl, sw , Slr );

					for ( int j = 0; j < l; j ++ )
					{
						Sl[l - j - 1] = Slr[j];
					}
				}
			}

			int * M = ( int * ) calloc ( mm, sizeof ( int ) );
			unsigned int a = 0;
			unsigned int b = 0;

			for ( int j = 0; j < mm; j++ )
			{
				if ( j < l )
					M[j] = Sl[a++];
				if ( j >= l + mf[jj] )
					M[j] = Sr[b++];
			}

                        int left  = max ( 0, ( int ) ( ind[jj] + mf[jj] - m[w] )  );
                        int right = min ( ind[jj], m[w] - 1 );

                        int dist = 0;

                        for ( int j = left; j <= right; j++ )
			{

				if ( j == left )
				{
					for ( int l = left; l < left + m[w]; l++ )
					{
						dist += M[l];
					}
				}
				else
				{
					dist = dist - M[j - 1] + M[j + m[w] - 1];
				}

				if ( dist <= sw . k )
				{
					int pos = ii - l + j;

					if ( ( * NOcc )[w] >= MOcc[w] )
					{
						( * POcc ) [w] = ( TPOcc * ) realloc ( ( * POcc )[w],   ( MOcc[w] + ALLOC_SIZE ) * sizeof ( TPOcc ) );
						MOcc[w] += ALLOC_SIZE;
					}
					( *POcc )[w][ ( * NOcc )[w] ] . pos = pos;	//text index
					( *POcc )[w][ ( * NOcc )[w] ] . err = dist;	
					( *POcc )[w][ ( * NOcc )[w] ] . rot = j;	
					( * NOcc )[w] = ( * NOcc )[w] + 1;
				}
			}	

			if ( r > 0 )
				free ( Sr );

			if ( l > 0 )
			{
				free ( Sl );
				free ( Slr );
			}

			free ( M );
		
			jj = d_occ[jj];
			
		} while ( jj != -1 );
	}

	free ( F );
	free ( P );

	/* A bit of post-processing to get rid of the duplicated occurrences */
	for ( i = 0; i < d; i++ )
	{
		vector<mytuple> data;

		int j;
		for ( j = 0; j < ( * NOcc )[i]; j++ )
  			data.push_back(make_tuple(( *POcc )[i][j] . pos, ( *POcc )[i][j] . err, ( *POcc )[i][j] . rot));

  		std::sort( data.begin(), data.end() );
		
		j = 0;
		for(vector<mytuple>::iterator iter = data.begin(); iter != data.end(); iter++)
		{
    			( *POcc )[i][j] . pos = get<0>(*iter);
			( *POcc )[i][j] . err = get<1>(*iter);
			( *POcc )[i][j] . rot = get<2>(*iter);
			j++;
  		}

		/* Here we remove the occurrences by imlpementing the C++ unique function */
		unsigned int ue = 0;
                TPOcc * Last =  unique ( &( ( *POcc )[i][0] ), &( ( *POcc )[i][( * NOcc )[i]] ) );
	
		TPOcc * Current;
		for ( Current = &( ( *POcc )[i][0] ); Current != Last; ++Current, ++ ue );
                ( * NOcc )[i] = ue;

		/* Here we resize the vector to save space */
		if ( ( ( * NOcc )[i] ) > 0  )
			( * POcc ) [i]   = ( TPOcc * ) realloc ( ( * POcc ) [i], ( ( * NOcc )[i] ) * sizeof ( TPOcc ) );
	} 

        free ( MOcc );

	free ( m );
	for ( int i = 0; i < d; i++ )
	{
		free( xx[i] );
		free( xxr[i] );
	}	
	free ( xx );
	free ( xxr );

	free ( tr );

	free ( mf );
	free ( ind );
	free ( whe );

	for ( int i = 0; i < f * d; i++ )
		free ( seqs[i] );
	free ( seqs );

	free ( l_occ );
	free ( d_occ );
	free ( dups );

	return 1;
}

unsigned int r_errors_vec	( unsigned int ** D, unsigned int m, unsigned int n, struct TSwitch sw, int * Sr )
{
	int i, j;

	for ( i = 1; i < m + 1; i ++ )
	{
		/* We find the minimum value of row i */
		int min = m;
		
		for ( j = 1; j < n + 1; j++ )
		{
			if ( D[i][j] < min )
				min = D[i][j];
		}

		/* We store the max value of row i */
		Sr[i - 1] = min;
	}

	return ( 1 );
}

unsigned int l_errors_vec	( unsigned int ** D, unsigned int m, unsigned int n, int ** H, struct TSwitch sw, int * Sl, int * Pl )
{
	int i, j, k;

	unsigned int i_max = min ( m, n + sw . k );

	/* We will process each row of the DP matrix*/
        for ( i = 1; i < m + 1; i ++ )
	{
		/* We find the minimum value of row i */
		int min = m;
		unsigned int jj = n;

	        for ( j = 1; j < n + 1; j++ )
		{
			if ( D[i][j] < min )
			{
				min = D[i][j];
				jj = j;
			}
		}

		/* We store the max value of row i */
		Sl[i - 1] = min;

		/* We count the number of gaps inserted IN THE TEXT ONLY for this value of the DP matrix: D[i,jj] */
		k = i;
		j = jj;
		int txt_chars = 0;
		while ( !( k == 0 && j == 0 ) )
		{
			if ( H[k][j] == 0 )
			{
				j = j - 1;	
				k = k - 1;
				//there is no gap
			}
			if ( H[k][j] == -1 )
			{
				k = k - 1;
				txt_chars--;
				//a gap is inserted in the text so we decrease the number of txt chars participating in the alignment
			}
			if ( H[k][j] == 1 )
			{
				j = j - 1;	
				txt_chars++;
				//a gap is inserted in the pattern so we increase the number of txt chars participating in the alignment
			}
		}

		/* We store the number of text chars participating in the alignment */
		Pl[i - 1] = txt_chars;
	}

	return ( 1 );
}

/* Edit distance */
unsigned int edit_distance ( unsigned char * x, unsigned int m, unsigned char * y, unsigned int n, struct TSwitch  sw, unsigned int ** D )
{
	unsigned int i, j;
	unsigned int sub = MIS_PEN;
	unsigned int gap = GAPS_PEN;

	for ( i = 1; i < m + 1; i++ )
    	{
      		D[i][0] = i * gap;
    	}
	for ( j = 1; j < n + 1; j++ )
    	{
      		D[0][j] = j * gap;
    	}

	for ( j = 1; j < n + 1; j++ )
    	{
		for ( i = 1; i < m + 1; i++ )
        	{
          		if ( x[i - 1] == y[j - 1] )
			{
            			D[i][j] = D[i - 1][j - 1];      //a match
			}
          		else
			{
				if ( D[i-1][j-1] + sub <= D[i][j-1] + gap && D[i-1][j-1] + sub <= D[i-1][j] + gap )
				{
					D[i][j] = D[i-1][j-1] + sub;
				}
				else if ( D[i-1][j] + gap <= D[i][j-1] + gap && D[i-1][j] + gap <= D[i-1][j-1] + sub )
				{
					D[i][j] = D[i - 1][j] + gap;
				}
				else if ( D[i][j-1] + gap <= D[i - 1][j] + gap && D[i][j-1] + gap <= D[i - 1][j-1] + sub )
				{
					D[i][j] = D[i][j-1] + gap;
				}
			}
        	}
    	}
	return ( 1 );
}  

/* Edit distance with backtracking */
unsigned int edit_distance_wbt ( unsigned char * x, unsigned int m, unsigned char * y, unsigned int n, struct TSwitch  sw, unsigned int ** D, int ** H )
{
	unsigned int i, j;
	unsigned int sub = MIS_PEN;
	unsigned int gap = GAPS_PEN;

	for ( i = 1; i < m + 1; i++ )
    	{
      		D[i][0] = i * gap;
      		H[i][0] = -1;
    	}
	for ( j = 1; j < n + 1; j++ )
    	{
      		D[0][j] = j * gap;
      		H[0][j] = 1;
    	}

	for ( j = 1; j < n + 1; j++ )
    	{
		for ( i = 1; i < m + 1; i++ )
        	{
          		if ( x[i - 1] == y[j - 1] )
			{
            			D[i][j] = D[i - 1][j - 1];      //a match
            			H[i][j] = 0;       	     	
			}
          		else
			{
		
				if ( D[i-1][j-1] + sub <= D[i][j-1] + gap && D[i-1][j-1] + sub <= D[i-1][j] + gap )
				{
					D[i][j] = D[i-1][j-1] + sub;
            				H[i][j] = 0;       	//a substitution
				}
				else if ( D[i-1][j] + gap <= D[i][j-1] + gap && D[i-1][j] + gap <= D[i-1][j-1] + sub )
				{
					D[i][j] = D[i - 1][j] + gap;
            				H[i][j] = -1;   	//a deletion
				}
				else if ( D[i][j-1] + gap <= D[i - 1][j] + gap && D[i][j-1] + gap <= D[i - 1][j-1] + sub )
				{
					D[i][j] = D[i][j - 1] + gap;
            				H[i][j] = 1; 		//an insertion
				}
			}
        	}
    	}
	return ( 1 );
}  

/* Hamming distance */
unsigned int hamming_distance ( unsigned char * x, unsigned int m, unsigned char * y, unsigned int n, struct TSwitch sw, unsigned int * S )
{
	for ( unsigned int i = 0; i < m; i++ )
    	{
		if ( i < n )
		{
			if ( x[i] != y[i] )
      				S[i] = 1;
			else
				S[i] = 0;
		}
		else
			S[i] = m - 1;	
    	}
	return ( 1 );
}  
