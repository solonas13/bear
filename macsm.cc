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
#include <vector>
#include <tuple>
#include <algorithm>

#include "beardefs.h"
#include "filter.h"
#include "aca.h"
#include "EDNAFULL.h"
#include "EBLOSUM62.h"
#include "globals.h"

#define OPT_LEN_FRA 20

using namespace std;

/* Definition of global variables defined in globals.h */
int * gF  = NULL;
int * gP  = NULL;
int gMatches = 0;
int gMax_alloc_matches = ALLOC_SIZE;

double gettime( void )
{
    struct timeval ttime;
    gettimeofday( &ttime , 0 );
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
}

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

			double * Sr;
			double * Sl;
			double * Slr;
			int * Pl;
			int * Plr;
			
			/* the right part*/
			if ( r > 0 )
			{
				Sr = ( double * ) calloc ( r,  sizeof ( double ) );

				double ** D;
				D = ( double ** ) calloc ( r + 1,  sizeof ( double * ) );
				for ( int j = 0; j < r + 1; j ++ )
					D[j] = ( double * ) calloc ( exr + 1,  sizeof ( double ) );

				nw_algorithm ( &xx[w][pxx], r, &t[pt], exr, sw, D );
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
				Sl = ( double * ) calloc ( l,  sizeof ( double ) );
				Pl = ( int * ) calloc ( l,  sizeof ( int ) );
				Slr = ( double* ) calloc ( l,  sizeof ( double ) );
				Plr = ( int * ) calloc ( l,  sizeof ( int ) );

				double ** D;
				int ** H;

				D = ( double ** ) calloc ( l + 1,  sizeof ( double * ) );
				for ( int j = 0; j < l + 1; j ++ )
					D[j] = ( double * ) calloc ( exl + 1,  sizeof ( double ) );

				H = ( int ** ) calloc ( l + 1,  sizeof ( int * ) );
				for ( int j = 0; j < l + 1; j ++ )
					H[j] = ( int * ) calloc ( exl + 1,  sizeof ( int ) );

				nw_algorithm_wbt ( &xxr[w][pxxr], l, &tr[ptr], exl, sw, D, H );
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

			double * M = ( double * ) calloc ( mm, sizeof ( double ) );
			unsigned int a = 0;
			unsigned int b = 0;

			for ( int j = 0; j < mm; j++ )
			{
				if ( j < l )
					M[j] = Sl[a++];
				if ( j >= l + mf[jj] )
					M[j] = Sr[b++];
			}		

			//double min_sim = ( double ) m[w] * MIN_SIM * MATCH;

			double frag_score = 0;
			for ( int j = 0; j < mf[jj]; j++ )
				frag_score += nuc_delta ( t[ii + j], t[ii + j] );

                        int left  = max ( 0, ( int ) ( ind[jj] + mf[jj] - m[w] )  );	//leftmost rotation j
                        int right = min ( ind[jj], m[w] - 1 );				//rightmost rotation j

			for ( int j = left; j <= right; j++ )
			{
				double score = 0;
				double led = 0;
				double red = 0;
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

				score = led + frag_score + red;

				if ( score >= sw . min_sim )
				{
					int txt_chars = l + pl;
					int pos = ii - txt_chars + j;

					if ( ( * NOcc )[w] >= MOcc[w] )
					{
						( * POcc ) [w] = ( TPOcc * ) realloc ( ( * POcc )[w],   ( MOcc[w] + ALLOC_SIZE ) * sizeof ( TPOcc ) );
						MOcc[w] += ALLOC_SIZE;
					}
					( *POcc )[w][ ( * NOcc )[w] ] . pos = pos;	//text index
					( *POcc )[w][ ( * NOcc )[w] ] . err = score;	
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

                        int mis = 0;

                        for ( int j = left; j <= right; j++ )
			{

				if ( j == left )
				{
					for ( int l = left; l < left + m[w]; l++ )
					{
						mis += M[l];
					}
				}
				else
				{
					mis = mis - M[j - 1] + M[j + m[w] - 1];
				}

				if ( mis <= sw . k )
				{
					int pos = ii - l + j;

					if ( ( * NOcc )[w] >= MOcc[w] )
					{
						( * POcc ) [w] = ( TPOcc * ) realloc ( ( * POcc )[w],   ( MOcc[w] + ALLOC_SIZE ) * sizeof ( TPOcc ) );
						MOcc[w] += ALLOC_SIZE;
					}
					( *POcc )[w][ ( * NOcc )[w] ] . pos = pos;	//text index
					( *POcc )[w][ ( * NOcc )[w] ] . err = mis;	
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

TPOcc * unique ( TPOcc * first, TPOcc * last )
{
  	if ( first == last ) return last;
	
  	TPOcc * result = first;
  	while (  ( ++first != last ) )
  	{
    		if ( ! ( ( * result ) . pos == ( * first ) . pos && ( * result ) . err == ( * first ) . err && ( * result ) . rot == ( * first ) . rot ) )
		{
			++result;
      			( * result ) . pos = ( * first ) . pos;
      			( * result ) . err = ( * first ) . err;
      			( * result ) . rot = ( * first ) . rot;
		}
  	}
	++result;
  	return ( result );
}

unsigned int create_rotation ( unsigned char * x, unsigned int offset, unsigned char * rotation )
{
	unsigned int m = strlen ( ( char * ) x );
	memmove ( &rotation[0], &x[offset], m - offset );
	memmove ( &rotation[m - offset], &x[0], offset );
	rotation[m] = '\0';
	return 1;
}

int binSearch( int value, TPat * array, int num )
{
   int mid = (num - 1) / 2;
   int index = -1;

   if ( num > 0 )
   {
      if ( value < array[mid] . start )
      {
         index = binSearch(value, array, mid);
      }
      else if ( value == array[mid] . start )
      {
         index = mid;
      }
      else if (value > array[mid] . start )
      {
         index = binSearch(value, array + mid + 1, num - mid - 1);

         if ( index >= 0 )
         {
            index += (mid + 1);
         }
      }
   }

   return index;
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

unsigned int r_errors_vec	( double ** D, unsigned int m, unsigned int n, struct TSwitch sw, double * Sr )
{
	int i, j;

	/* We will process each row of the DP matrix*/
	unsigned int i_max = min ( m, n + sw . gaps );

//	for ( i = 1; i < m + 1; i ++ )
	for( i = 1; i < i_max + 1; i++)
	{
		/* We find the minimum value of row i */
		double max = -DBL_MAX;
		
		//for ( j = 1; j < n + 1; j++ )
		unsigned int j_min = max ( 1, ( int ) ( i - sw . gaps ));
		unsigned int j_max = min ( n, ( int ) ( i + sw . gaps ));
		//fprintf( stderr, "%d %d %d %d %d %d", sw . gaps, j_min, j_max, n, i, m ); getchar();	
		for( j = j_min; j <= j_max; j++ )
		{
			if ( D[i][j] > max )
				max = D[i][j];
		}
		//fprintf( stderr, "Ok!" ); getchar();	

		/* We store the max value of row i */
		Sr[i - 1] = max;
	}

	return ( 1 );
}

unsigned int l_errors_vec	( double ** D, unsigned int m, unsigned int n, int ** H, struct TSwitch sw, double * Sl, int * Pl )
{
	int i, j, k;

	unsigned int i_max = min ( m, n + sw . gaps );

	/* We will process each row of the DP matrix*/
//        for ( i = 1; i < m + 1; i ++ )
	for( i = 1; i < i_max + 1; i++)
	{
		/* We find the minimum value of row i */
		double max = -DBL_MAX;
		unsigned int jj = n;

		unsigned int j_min = max ( 1, ( int ) ( i - sw . gaps ));
		unsigned int j_max = min ( n, ( int ) ( i + sw . gaps ));
		
//	        for ( j = 1; j < n + 1; j++ )
		for( j = j_min; j <= j_max; j++ )
		{
			if ( D[i][j] > max )
			{
				max = D[i][j];
				jj = j;
			}
		}

		/* We store the max value of row i */
		Sl[i - 1] = max;

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

unsigned int nw_algorithm ( unsigned char * p, unsigned int m, unsigned char * t, unsigned int n, struct TSwitch sw, double ** T )
{
	int i;
	int j;
        double **       D;
        double **       I;
	double g = sw . O;
        double h = sw . E;
	double matching_score = 0;
	unsigned int j_max = min ( n, m + sw . gaps );

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

        for( j = 1; j < j_max + 1; j++)
        {
                unsigned int i_min = max ( 1, ( int ) ( j - sw . gaps ));
                unsigned int i_max = min ( m, ( int ) ( j + sw . gaps ));
                for( i = i_min; i <= i_max; i++ )
                {

			D[i][j] = max ( D[i - 1][j] + h, T[i - 1][j] + g );
			double u = D[i][j];
			I[i][j] = max ( I[i][j - 1] + h, T[i][j - 1] + g );
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
	unsigned int j_max = min ( n, m + sw . gaps );

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

        for( j = 1; j < j_max + 1; j++)
        {
                unsigned int i_min = max ( 1, ( int ) ( j - sw . gaps ));
                unsigned int i_max = min ( m, ( int ) ( j + sw . gaps ));
                for( i = i_min; i <= i_max; i++ )
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

/* Returns the score for matching character a and b based on EDNAFULL matrix */
int nuc_delta ( char a, char b )
 {
   unsigned int index_a = nuc_char_to_index ( a );
   unsigned int index_b = nuc_char_to_index ( b );

   if ( ( index_a < NUC_SCORING_MATRIX_SIZE ) && ( index_b < NUC_SCORING_MATRIX_SIZE ) )
     return ( EDNAFULL_matrix[ index_a ][ index_b ] );
   else //Error
     return ( ERR );
 }

/* Returns the index of char a in EDNAFULL matrix */
unsigned int nuc_char_to_index ( char a )
 {
   unsigned int index; 

   switch ( a )
    {
      case 'A':
        index = 0; break;

      case 'T':
        index = 1; break;

      case 'G':
        index = 2; break;

      case 'C':
        index = 3; break;

      case 'S':
        index = 4; break;

      case 'W':
        index = 5; break;

      case 'R':
        index = 6; break;

      case 'Y':
        index = 7; break;

      case 'K':
        index = 8; break;

      case 'M':
        index = 9; break;

      case 'B':
        index = 10; break;

      case 'V':
        index = 11; break;

      case 'H':
        index = 12; break;

      case 'D':
        index = 13; break;

      case 'N':
        index = 14; break;

      default:
        fprintf ( stderr, "Error: unrecognizable character in one of the nucleotide sequences!!!\n" );
        index = ERR; break;
    }
   
   return ( index );
 }

/* Returns the index of char a in EBLOSUM62 matrix */
unsigned int pro_char_to_index ( char a )
{
	unsigned int index;
	switch ( a )
	{
		case 'A':
		index = 0; break;
		case 'R':
		index = 1; break;
		case 'N':
		index = 2; break;
		case 'D':
		index = 3; break;
		case 'C':
		index = 4; break;
		case 'Q':
		index = 5; break;
		case 'E':
		index = 6; break;
		case 'G':
		index = 7; break;
		case 'H':
		index = 8; break;
		case 'I':
		index = 9; break;
		case 'L':
		index = 10; break;
		case 'K':
		index = 11; break;
		case 'M':
		index = 12; break;
		case 'F':
		index = 13; break;
		case 'P':
		index = 14; break;
		case 'S':
		index = 15; break;
		case 'T':
		index = 16; break;
		case 'W':
		index = 17; break;
		case 'Y':
		index = 18; break;
		case 'V':
		index = 19; break;
		case 'B':
		index = 20; break;
		case 'Z':
		index = 21; break;
		case 'X':
		index = 22; break;
		case '*':
		index = 23; break;
		default:
		fprintf ( stderr, "Error: unrecognizable character in one of the protein sequences!!!\n" );
		index = ERR; break;
	}
	return ( index );
}
/* Returns the score for matching character a and b based on EBLOSUM62 matrix */
int pro_delta ( char a, char b )
{
	unsigned int index_a = pro_char_to_index( a );
	unsigned int index_b = pro_char_to_index( b );
	if ( ( index_a < PRO_SCORING_MATRIX_SIZE ) && ( index_b < PRO_SCORING_MATRIX_SIZE ) )
		return ( EBLOSUM62_matrix[ index_a ][ index_b ] );
	else //Error
	return ( ERR );
}
