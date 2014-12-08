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

#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <float.h>
#include <limits.h>
#include <omp.h>
#include "beardefs.h"

int main(int argc, char **argv)
{

	struct TSwitch  sw;

	FILE *          in_fd;                  // the input file descriptor
	FILE *          out_fd;                 // the output file descriptor
	FILE *          outl_fd;                // the outliers file descriptor
        char *          patterns_filename;      // the patterns input file name
        char *          text_filename;          // the text input file name
        char *          output_filename;        // the output file name
        char *          rotations_filename;     // the rotations file name
        char *          outliers_filename;      // the outliers file name
	struct TPat   * pat     = NULL;
        unsigned char ** seq    = NULL;         // the sequence(s) in memory
        unsigned char ** seq_id = NULL;         // the sequence(s) id in memory
        unsigned char *  t_id   = NULL;         // the text id in memory
        unsigned char *  t      = NULL;         // the text in memory
	char *          alphabet;               // the alphabet
	unsigned int    i, j;
	unsigned int    d;	
	unsigned int    total_length = 0;

	/* Decodes the arguments */
        i = decode_switches ( argc, argv, &sw );

	/* Check the arguments */
        if ( i < 7 )
        {
                usage ();
                return ( 1 );
        }
        else
        {
                if      ( ! strcmp ( "DNA", sw . alphabet ) )   { alphabet = ( char * ) DNA;  sw . matrix = 0; }
                else if ( ! strcmp ( "PROT", sw . alphabet ) )  { alphabet = ( char * ) PROT; sw . matrix = 1; }
                else
                {
                        fprintf ( stderr, " Error: alphabet argument a should be `DNA' for nucleotide sequences or `PROT' for protein sequences or `USR' for sequences over a user-defined alphabet!\n" );
                        return ( 1 );
                }

		d       = sw . d;
		if ( d > 2 )	
		{
                	usage ();
			return ( 1 );
		}

		if ( d == 0 )
		{
			if ( sw . D > 1 )
			{
				fprintf ( stderr, " Error: The distance model (-D) must be either 0 or 1!\n" );
				return ( 1 );
			}

		}

		if ( d == 1 )
		{
			if ( sw . w > WORD_LEN - 1 )
			{
				fprintf ( stderr, " Error: The length of factor (-w) must be smaller than %d!\n", WORD_LEN );
				return ( 1 );
			}
			if ( sw . w <= 0 )
			{
				fprintf ( stderr, " Error: The length of factor (-w) must be larger than 1!\n" );
				return ( 1 );
			}
			if ( sw . D > 1 )
			{
				fprintf ( stderr, " Error: The distance model (-D) must be either 0 or 1!\n" );
				return ( 1 );
			}

			if ( sw . k >= sw . w )
			{
				fprintf ( stderr, " Error: The maximum distance (-k) must be smaller than the factor length (-w)!\n" );
				return ( 1 );
			}
		}

		if ( d == 2 )
		{
			if ( sw . R <= 0 || sw . R >= 1 )
			{
				fprintf ( stderr, " Error: The similarity ratio (-R) must be in (0,1)!\n" );
				return ( 1 );
			}
			if ( sw . O >= 0 )
			{
				fprintf ( stderr, " Error: The gap opening penalty (-O) must be < 0!\n" );
				return ( 1 );
			}
			if ( sw . E >= 0 )
			{
				fprintf ( stderr, " Error: The gap extension penalty (-E) must be < 0!\n" );
				return ( 1 );
			}
			if ( sw . A > 1 )
			{
				fprintf ( stderr, " Error: The alignment type (-A) must be either 0 or 1!\n" );
				return ( 1 );
			}
		}


                patterns_filename       = sw . patterns_filename;
		if ( patterns_filename == NULL )
		{
			fprintf ( stderr, " Error: Cannot open file for input!\n" );
			return ( 1 );
		}
		output_filename         = sw . output_filename;
		
                text_filename           = sw . text_filename;
                rotations_filename      = sw . rotations_filename;
                outliers_filename       = sw . outliers_filename;

                if ( text_filename != NULL && rotations_filename != NULL  )
		{
			fprintf ( stderr, " Error: Text filename cannot be used with rotations filename!\n" );
			return ( 1 );
		}
                if ( outliers_filename != NULL && rotations_filename != NULL  )
		{
			fprintf ( stderr, " Error: Outliers filename cannot be used with rotations filename!\n" );
			return ( 1 );
		}
                if ( outliers_filename != NULL && text_filename != NULL  )
		{
			fprintf ( stderr, " Error: Outliers filename cannot be used with text filename!\n" );
			return ( 1 );
		}

		omp_set_num_threads( sw . T );
        }


	double start = gettime();

	/* Read the (Multi)FASTA file in memory */
	fprintf ( stderr, " Reading the (Multi)FASTA input file: %s\n", patterns_filename );
	if ( ! ( in_fd = fopen ( patterns_filename, "r") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", patterns_filename );
		return ( 1 );
	}

	char c;
        unsigned int num_seqs = 0;           // the total number of sequences considered
	unsigned int max_alloc_seq_id = 0;
	unsigned int max_alloc_seq = 0;
	c = fgetc( in_fd );
	do
	{
		if ( c != '>' )
		{
			fprintf ( stderr, " Error: input file %s is not in FASTA format!\n", patterns_filename );
			return ( 1 );
		}
		else
		{
			if ( num_seqs >= max_alloc_seq_id )
			{
				seq_id = ( unsigned char ** ) realloc ( seq_id,   ( max_alloc_seq_id + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
				max_alloc_seq_id += ALLOC_SIZE;
			}

			unsigned int max_alloc_seq_id_len = 0;
			unsigned int seq_id_len = 0;

			seq_id[ num_seqs ] = NULL;

			while ( ( c = fgetc( in_fd ) ) != EOF && c != '\n' )
			{
				if ( seq_id_len >= max_alloc_seq_id_len )
				{
					seq_id[ num_seqs ] = ( unsigned char * ) realloc ( seq_id[ num_seqs ],   ( max_alloc_seq_id_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
					max_alloc_seq_id_len += ALLOC_SIZE;
				}
				seq_id[ num_seqs ][ seq_id_len++ ] = c;
			}
			seq_id[ num_seqs ][ seq_id_len ] = '\0';
			
		}

		if ( num_seqs >= max_alloc_seq )
		{
			seq = ( unsigned char ** ) realloc ( seq,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
			max_alloc_seq += ALLOC_SIZE;
		}

		unsigned int seq_len = 0;
		unsigned int max_alloc_seq_len = 0;

		seq[ num_seqs ] = NULL;

		while ( ( c = fgetc( in_fd ) ) != EOF && c != '>' )
		{
			if( seq_len == 0 && c == '\n' )
			{
				fprintf ( stderr, " Omitting empty sequence in file %s!\n", patterns_filename );
				c = fgetc( in_fd );
				break;
			}
			if( c == '\n' || c == ' ' ) continue;

			c = toupper( c );

			if ( seq_len >= max_alloc_seq_len )
			{
				seq[ num_seqs ] = ( unsigned char * ) realloc ( seq[ num_seqs ],   ( max_alloc_seq_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
				max_alloc_seq_len += ALLOC_SIZE;
			}

			if( strchr ( alphabet, c ) )
			{
				seq[ num_seqs ][ seq_len++ ] = c;
			}
			else
			{
				fprintf ( stderr, " Error: input file %s contains an unexpected character %c!\n", patterns_filename, c );
				return ( 1 );
			}

		}

		if( seq_len != 0 )
		{
			if ( seq_len >= max_alloc_seq_len )
			{
				seq[ num_seqs ] = ( unsigned char * ) realloc ( seq[ num_seqs ],   ( max_alloc_seq_len + ALLOC_SIZE ) * sizeof ( unsigned char ) ); 
				max_alloc_seq_len += ALLOC_SIZE;
			}
			seq[ num_seqs ][ seq_len ] = '\0';
			total_length += seq_len;
			num_seqs++;
		}
		
	} while( c != EOF );

	seq[ num_seqs ] = NULL;

	pat = ( TPat * ) calloc ( num_seqs , sizeof ( TPat ) ); 
	if( ( pat == NULL) ) 
	{
		fprintf( stderr, " Error: Cannot allocate memory!\n" );
		exit( EXIT_FAILURE );
	}

	for ( i = 0; i < num_seqs; i ++ )
	{
		if ( i == 0 )	pat[i] . start = 0;
		else 		pat[i] . start  += pat[i - 1] . start + strlen ( ( char * ) seq[ i - 1 ] );
		pat[i] . offset = 0;
	}

	if ( fclose ( in_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	if ( rotations_filename != NULL )
	{
		if ( ! ( out_fd = fopen ( rotations_filename, "w") ) )
		{
			fprintf ( stderr, " Error: Cannot open file %s!\n", rotations_filename );
			return ( 1 );
		}
		for ( i = 0; i < num_seqs; i ++ )
		{
			unsigned char * tmp;
			unsigned int l = strlen ( ( char * ) seq[i] );
			tmp = ( unsigned char * ) malloc ( ( l + 1 ) * sizeof ( unsigned char ) );
			if( ( tmp == NULL) ) 
			{
				fprintf( stderr, " Error: Cannot allocate memory!\n" );
				exit( EXIT_FAILURE );
			}
			int cp = rand() % ( l ) + 1;
			memcpy ( tmp, &seq[i][cp], l - cp );
			memcpy ( &tmp[l - cp], &seq[i][0], cp ); 
			tmp[l] = '\0';
			fprintf( out_fd, ">%s\n", seq_id[i]);
			fprintf( out_fd, "%s\n", tmp);
			free ( tmp );
		}
		if ( fclose ( out_fd ) )
		{
			fprintf( stderr, " Error: file close error!\n");
			return ( 1 );
		}
	}
	else if ( text_filename != NULL )
	{
		fprintf ( stderr, " Reading the FASTA input file: %s\n", text_filename );
		if ( ! ( in_fd = fopen ( text_filename, "r") ) )
		{
			fprintf ( stderr, " Error: Cannot open file %s!\n", text_filename );
			return ( 1 );
		}

		c = fgetc( in_fd );
		if ( c != '>' )
		{
			fprintf ( stderr, " Error: input file %s is not in FASTA format!\n", text_filename );
			return ( 1 );
		}
		else
		{
			unsigned int max_alloc_t_id_len = 0;
			unsigned int t_id_len = 0;

			while ( ( c = fgetc( in_fd ) ) != EOF && c != '\n' )
			{
				if ( t_id_len >= max_alloc_t_id_len )
				{
					t_id = ( unsigned char * ) realloc ( t_id,   ( max_alloc_t_id_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
					max_alloc_t_id_len += ALLOC_SIZE;
				}
				t_id[ t_id_len++ ] = c;
			}
			t_id[ t_id_len ] = '\0';
			
		}

		unsigned int t_len = 0;
		unsigned int max_alloc_t_len = 0;

		while ( ( c = fgetc( in_fd ) ) != EOF )
		{
			if( t_len == 0 && c == '\n' )
			{
				fprintf ( stderr, " Error: empty sequence in file %s!\n", text_filename );
				return ( 1 );
			}
			if( c == '\n' ) continue;

			c = toupper( c );

			if ( t_len >= max_alloc_t_len )
			{
				t = ( unsigned char * ) realloc ( t,   ( max_alloc_t_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
				max_alloc_t_len += ALLOC_SIZE;
			}

			if( strchr ( alphabet, c ) )
			{
				t[ t_len++ ] = c;
			}
			else
			{
				fprintf ( stderr, " Error: input file %s contains an unexpected character %c!\n", text_filename, c );
				return ( 1 );
			}

		}

		if( t_len != 0 )
		{
			if ( t_len >= max_alloc_t_len )
			{
				t = ( unsigned char * ) realloc ( t,   ( max_alloc_t_len + ALLOC_SIZE ) * sizeof ( unsigned char ) ); 
				max_alloc_t_len += ALLOC_SIZE;
			}
			t[ t_len ] = '\0';
		}

		if ( fclose ( in_fd ) )
		{
			fprintf( stderr, " Error: file close error!\n");
			return ( 1 );
		}

		/* Allocate the occurences structure */
		TPOcc ** POcc = NULL;
		POcc = ( TPOcc ** ) calloc ( num_seqs , sizeof ( TPOcc * ) ); 
		if( ( POcc == NULL) ) 
		{
			fprintf( stderr, " Error: Cannot allocate memory!\n" );
			exit( EXIT_FAILURE );
		}

		for ( i = 0; i < num_seqs; i ++ )
			POcc[i] = NULL;

		unsigned int * NOcc;
		NOcc = ( unsigned int * ) calloc ( num_seqs , sizeof ( unsigned int ) );
		if( ( NOcc == NULL) )
		{
			fprintf( stderr, " Error: Cannot allocate memory!\n" );
			exit( EXIT_FAILURE );
		}
		
		/* Multiple Circular Approximate String Matching */
		fprintf ( stderr, " Starting the multiple circular approximate string matching\n" );
		if ( d == 0 )
		{
			if ( ! ( macsmf_hd ( seq, t, sw, &POcc, &NOcc ) ) )
			{
				fprintf( stderr, " Error: macsmf_ms() failed!\n" );
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			if ( ! ( macsmf_ed ( seq, t, sw, &POcc, &NOcc ) ) )
			{
				fprintf( stderr, " Error: macsmf_ms() failed!\n" );
				exit(EXIT_FAILURE);
			}
		}

		fprintf ( stderr, " Preparing the output\n" );
		if ( ! ( out_fd = fopen ( output_filename, "w") ) )
		{
			fprintf ( stderr, " Error: Cannot open file %s!\n", output_filename );
			return ( 1 );
		}

		for ( i = 0; i < num_seqs; i ++ )
		{
			for ( int j = 0; j < NOcc[i]; j ++ )
			{
				if ( d == 0 )
					fprintf( out_fd, "%s %s %d %d %d\n", seq_id[i], seq[i], POcc[i][j] . pos, ( int ) POcc[i][j] . err, POcc[i][j] . rot );
				else
					fprintf( out_fd, "%s %s %d %lf %d\n", seq_id[i], seq[i], POcc[i][j] . pos, POcc[i][j] . err, POcc[i][j] . rot );
			}
		}

		
		if ( fclose ( out_fd ) )
		{
			fprintf( stderr, " Error: file close error!\n");
			return ( 1 );
		}

		free ( t_id );
		free ( t );
		free ( NOcc );
		for ( i = 0; i < num_seqs; i ++ )
			free ( POcc[i] );
		free ( POcc);
	}
	else
	{
		TPOcc ** D;

		int * R;
		R = ( int * ) calloc ( num_seqs , sizeof ( int ) );

		if ( d == 0 ) 
		{
			fprintf ( stderr, " Comparison model: small pairwise distance.\n" );
			if ( ( D = ( TPOcc ** ) calloc ( ( num_seqs ) , sizeof( TPOcc * ) ) ) == NULL )
			{
				fprintf( stderr, " Error: Cannot allocate memory!\n" );
				exit( EXIT_FAILURE );
			}
			for ( int i = 0; i < num_seqs; i++ )
			{
				if ( ( D[i] = ( TPOcc * ) calloc ( ( num_seqs ) , sizeof( TPOcc ) ) ) == NULL )
				{
					fprintf( stderr, " Error: Cannot allocate memory!\n" );
					exit( EXIT_FAILURE );
				}
			}
			if ( sw . D == 0 )
			{
				fprintf ( stderr, " Distance model: Hamming distance with k = %d.\n", sw . k );
				//#pragma omp parallel for
				for ( int i = 0; i < num_seqs; i++ )
				{
					unsigned int m = strlen ( ( char * ) seq[i] );
					if ( ( double ) sw . k / m > 0.10 )
						fprintf ( stderr, " Warning: Error ratio  %lf is too high for this method.\n", ( double ) sw . k / m );
					for ( int j = 0; j < num_seqs; j ++ )
					{
						if ( i == j ) continue;

						unsigned int n = strlen ( ( char * ) seq[j] );

						/* Initialise the arrays */
						D[i][j] . err = DBL_MAX;
						unsigned int distance = DBL_MAX;
						unsigned int rotation = 0;

						pcsa_hd ( seq[i], seq[j], sw, &rotation, &distance );
						
						D[i][j] . err = distance;
						D[i][j] . rot = rotation;		
					}
				}
			}
			else
			{
				fprintf ( stderr, " Distance model: edit distance with k = %d.\n", sw . k );
				//#pragma omp parallel for
				for ( int i = 0; i < num_seqs; i++ )
				{
					unsigned int m = strlen ( ( char * ) seq[i] );
					if ( ( double ) sw . k / m > 0.10 )
						fprintf ( stderr, " Warning: Error ratio  %lf is too high for this method.\n", ( double ) sw . k / m );
					for ( int j = 0; j < num_seqs; j ++ )
					{
						if ( i == j ) continue;

						unsigned int n = strlen ( ( char * ) seq[j] );

						/* Initialise the arrays */
						D[i][j] . err = DBL_MAX;
						unsigned int distance = DBL_MAX;
						unsigned int rotation = 0;

						pcsa_ed ( seq[i], seq[j], sw, &rotation, &distance );
						
						D[i][j] . err = distance;
						D[i][j] . rot = rotation;		
					}
				}
			}

			/* For every sequence i */

			#if 0
			/* Create the text */
			t = ( unsigned char * ) malloc ( ( total_length + 1 ) * sizeof ( unsigned char ) );
			unsigned int j = 0;
			for ( i = 0; i < num_seqs; i ++ )
			{
				unsigned int l = strlen ( ( char * ) seq[i] );
				memcpy ( &t[j], &seq[i][0], l );
				j += l;
			}
			t[ total_length ] = '\0';

			/* Allocate the occurences structure */
			TPOcc ** POcc = NULL;
			POcc = ( TPOcc ** ) calloc ( num_seqs , sizeof ( TPOcc * ) ); 
			if( ( POcc == NULL) ) 
			{
				fprintf( stderr, " Error: Cannot allocate memory!\n" );
				exit( EXIT_FAILURE );
			}

			for ( i = 0; i < num_seqs; i ++ )
				POcc[i] == NULL;

			unsigned int * NOcc;
			NOcc = ( unsigned int * ) calloc ( num_seqs , sizeof ( unsigned int ) );
			if( ( NOcc == NULL) )
			{
				fprintf( stderr, " Error: Cannot allocate memory!\n" );
				exit( EXIT_FAILURE );
			}
			
			/* Multiple Circular Approximate String Matching */
			if ( sw . D == 0 )
			{
				fprintf ( stderr, " Distance model: Hamming distance with k = %d.\n", sw . k );
				if ( ! ( macsmf_hd ( seq, t, sw, &POcc, &NOcc ) ) )
				{
					fprintf( stderr, " Error: macsmf_ms() failed!\n" );
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				fprintf ( stderr, " Distance model: edit distance with k = %d.\n", sw . k );
				if ( ! ( macsmf_ed ( seq, t, sw, &POcc, &NOcc ) ) )
				{
					fprintf( stderr, " Error: macsmf_ms() failed!\n" );
					exit(EXIT_FAILURE);
				}
			}

			D = ( TPOcc ** ) calloc ( num_seqs , sizeof ( TPOcc * ) );
			if( ( D == NULL) )
			{
				fprintf( stderr, " Error: Cannot allocate memory!\n" );
				exit( EXIT_FAILURE );
			}

			for ( i = 0; i < num_seqs; i ++ )
			{
				D[i] = ( TPOcc * ) calloc ( num_seqs , sizeof ( TPOcc ) );
				if( ( D[i] == NULL) )
				{
					fprintf( stderr, " Error: Cannot allocate memory!\n" );
					exit( EXIT_FAILURE );
				}

				/* Initialise the array */
				for ( j = 0; j < num_seqs; j++ )
					D[i][j] . err = INT_MAX;
			}

			for ( i = 0; i < num_seqs; i ++ )
			{
				for ( j = 0; j < NOcc[i]; j ++ )
				{
					int l = binSearch( POcc[i][j] . pos, pat, num_seqs );
					if ( l >= 0 && D[i][l] . err > POcc[i][j] . err )	
					{
						D[i][l] . err = POcc[i][j] . err;
						D[i][l] . pos = POcc[i][j] . pos;
						D[i][l] . rot = POcc[i][j] . rot;
					}
				}
			}

			free ( t );
			free ( NOcc );
			for ( i = 0; i < num_seqs; i ++ )
			{
				free ( POcc[i] );
				POcc[i] = NULL;
			}
			free ( POcc );
			#endif
		}

		if ( d == 1 )
		{	
			fprintf ( stderr, " Comparison model: pairwise distance of fixed-length factors.\n" );
			if ( ( D = ( TPOcc ** ) calloc ( ( num_seqs ) , sizeof( TPOcc * ) ) ) == NULL )
			{
				fprintf( stderr, " Error: Cannot allocate memory!\n" );
				exit( EXIT_FAILURE );
			}
			if ( sw . D == 0 )
			{
				fprintf ( stderr, " Distance model: Hamming distance with w = %d and k = %d.\n", sw . w, sw . k );
			}
			else
			{
				fprintf ( stderr, " Distance model: edit distance with w = %d and k = %d.\n", sw . w, sw . k );
			}

			/* For every sequence i */
			#pragma omp parallel for
			for ( int i = 0; i < num_seqs; i++ )
			{
				unsigned int m = strlen ( ( char * ) seq[i] );
				if ( ( D[i] = ( TPOcc * ) calloc ( ( num_seqs ) , sizeof( TPOcc ) ) ) == NULL )
				{
					fprintf( stderr, " Error: Cannot allocate memory!\n" );
					exit( EXIT_FAILURE );
				}

				for ( int ii = 0; ii < num_seqs; ii ++ )
				{
					if ( i == ii ) continue;

					TPOcc * M;

					/* Each row M[ii] stores the coordinates M . rot of the best match of length sw . w with every other sequence ii */
					if ( ( M = ( TPOcc * ) calloc ( ( m ) , sizeof( TPOcc ) ) ) == NULL )
					{
						fprintf( stderr, " Error: M could not be allocated!\n");
						exit ( 1 );
					}

					/* Here we compute these coordinates using the Max-Shift algorithm and store it to M[ii] */
					unsigned int n = strlen ( ( char * ) seq[ii] );
					unsigned char * xjw;
					xjw = ( unsigned char * ) calloc (  n + sw . w, sizeof( unsigned char ) );
					memmove( &xjw[0], seq[ii], n );
					memmove( &xjw[n], seq[ii], sw . w - 1 );
					xjw[n + sw . w - 1] = '\0';
					unsigned int nn = n + sw . w - 1;

					/* Initialise the arrays */
					for ( int jj = 0; jj < num_seqs; jj++ )
						M[jj] . err = m - 1;
					D[i][ii] . err = m  - 1;

					if ( sw . D == 0 )
						bcf_maxshift_hd_ls ( seq[i], m, xjw, nn, sw, M );
					else
						bcf_maxshift_ed_ls ( seq[i], m, xjw, nn, sw, M );

					unsigned int min_dist = sw . w;

					for ( int jj = sw . w - 1; jj < m; jj ++ )
					{
						if ( min_dist > M[jj] . err )
						{
							int iii = jj % m;
							int jjj = M[jj] . rot % n;

							min_dist = M[jj] . err;
							D[i][ii] . err = M[jj] . err;

							if ( iii >= jjj )
							{
								D[i][ii] . rot = iii - jjj;
							}
							else
							{
								int a = jjj - iii;
								int b = m - iii - 1;
								int c = min ( a, b );
								D[i][ii] . rot = m - c;		
							}
						}
					}
					free ( xjw );
					free ( M );
				}
			}	
		}

		if ( d == 2 )
		{	
			fprintf ( stderr, " Comparison model: affine gap penalty with sub matrices.\n" );
			if ( ( D = ( TPOcc ** ) calloc ( ( num_seqs ) , sizeof( TPOcc * ) ) ) == NULL )
			{
				fprintf( stderr, " Error: Cannot allocate memory!\n" );
				exit( EXIT_FAILURE );
			}

			/* Estimate the minimum pairwise allowed similarity */
			double total_score = 0;
			for ( i = 0; i < num_seqs; i ++ )
			{
				unsigned int l = strlen ( ( char * ) seq[i] );
				for ( int j = 0; j < l; j++ )
					total_score += ( sw . matrix ? pro_delta( seq[i][j], seq[i][j] ) : nuc_delta( seq[i][j], seq[i][j] ) ) ;
			}
			sw . min_sim = ( total_score * sw . R ) / ( num_seqs );

			if ( sw . A == 1 )
			{
				fprintf ( stderr, " Alignment type: local with O = %lf and E = %lf.\n", sw . O, sw . E );
				fprintf( stderr, " The minimum allowed similarity score is %lf.\n", sw . min_sim );

				/* For every sequence i */
				#pragma omp parallel for
				for ( int i = 0; i < num_seqs; i++ )
				{
					unsigned int m = strlen ( ( char * ) seq[i] );
					if ( ( D[i] = ( TPOcc * ) calloc ( ( num_seqs ) , sizeof( TPOcc ) ) ) == NULL )
					{
						fprintf( stderr, " Error: Cannot allocate memory!\n" );
						exit( EXIT_FAILURE );
					}

					for ( int j = 0; j < num_seqs; j ++ )
					{
						if ( i == j ) continue;

						TPOcc * M;

						/* Each row M[ii] stores the coordinates M . rot of the best match of length sw . w with every other sequence j */
						if ( ( M = ( TPOcc * ) calloc ( ( m ) , sizeof( TPOcc ) ) ) == NULL )
						{
							fprintf( stderr, " Error: M could not be allocated!\n");
							exit ( 1 );
						}

						/* Here we compute these coordinates using the Max-Shift algorithm and store it to M[ii] */
						unsigned int n = strlen ( ( char * ) seq[j] );
						unsigned char * xjw;
						xjw = ( unsigned char * ) calloc (  n + n, sizeof( unsigned char ) );
						memmove( &xjw[0], seq[j], n );
						memmove( &xjw[n], seq[j], n - 1 );
						xjw[n + n - 1] = '\0';
						unsigned int nn = n + n - 1;

						/* Initialise the arrays */
						for ( int jj = 0; jj < num_seqs; jj++ )
							M[jj] . err = -DBL_MAX;

						D[i][j] . err = -DBL_MAX;

						sw_ls ( seq[i], m, xjw, nn, sw, M );

						double max_sim = -DBL_MAX;

						for ( int jj = 0; jj < m; jj ++ )
						{
							if ( max_sim < M[jj] . err )
							{
								int iii = jj % m;
								int jjj = M[jj] . rot % n;

								max_sim = M[jj] . err;
								D[i][j] . err = M[jj] . err;

								if ( iii >= jjj )
								{
									D[i][j] . rot = iii - jjj;
								}
								else
								{
									int a = jjj - iii;
									int b = m - iii - 1;
									int c = min ( a, b );
									D[i][j] . rot = m - c;		
								}
							}
						}
						free ( xjw );
						free ( M );
					}
				}
			}
			else
			{	
				fprintf ( stderr, " Alignment type: global with O = %lf and E = %lf.\n", sw . O, sw . E );
				fprintf( stderr, " The minimum allowed similarity score is %lf.\n", sw . min_sim );

				for ( int i = 0; i < num_seqs; i++ )
				{
					unsigned int m = strlen ( ( char * ) seq[i] );
					if ( ( D[i] = ( TPOcc * ) calloc ( ( num_seqs ) , sizeof( TPOcc ) ) ) == NULL )
					{
						fprintf( stderr, " Error: Cannot allocate memory!\n" );
						exit( EXIT_FAILURE );
					}
				}

				/* For every sequence i */
				#pragma omp parallel for
				for ( int i = 0; i < num_seqs; i++ )
				{
					unsigned int m = strlen ( ( char * ) seq[i] );
					for ( int j = 0; j < num_seqs; j ++ )
					{
						if ( i == j ) continue;

						unsigned int n = strlen ( ( char * ) seq[j] );

						/* Initialise the arrays */
						D[i][j] . err = -DBL_MAX;
						double score = -DBL_MAX;
						int rot = 0;

						cyc_nw_ls ( seq[i], m, seq[j], n, sw, &score, &rot );
						
						D[i][j] . err = score;
						D[j][i] . rot = rot;		
					}
				}
			}	
		}

		#if 0
		for ( int i = 0; i < num_seqs; i ++ )
		{
			for ( int j = 0; j < num_seqs; j ++ )
			{
				fprintf( stderr, "%lf ", D[i][j] . err );
			}
			fprintf( stderr, "\n"  );
		}

		for ( int i = 0; i < num_seqs; i ++ )
		{
			for ( int j = 0; j < num_seqs; j ++ )
			{
				fprintf( stderr, "%d ", D[i][j] . rot );
			}
			fprintf( stderr, "\n"  );
		}
		#endif

		fprintf ( stderr, " Starting the clustering\n" );
		if ( d == 0 || d == 1 )
			upgma_dist ( D, num_seqs, sw, R, seq );
		if ( d == 2 )
			upgma_sim ( D, num_seqs, sw, R, seq );
		
		fprintf ( stderr, " Preparing the output\n" );

		if ( ! ( out_fd = fopen ( output_filename, "w") ) )
		{
			fprintf ( stderr, " Error: Cannot open file %s!\n", output_filename );
			return ( 1 );
		}

		unsigned int succ = 0;
		unsigned int unsucc = 0;
		for ( i = 0; i < num_seqs; i ++ )
		{
			if ( R[i] >= 0 )
			{
				R[i] = R[i] % strlen ( ( char * ) seq[i] );
				unsigned char * rotation;
				unsigned int m = strlen ( ( char * ) seq[i] );
				rotation = ( unsigned char * ) calloc ( m + 1, sizeof ( unsigned char ) );
				create_rotation ( seq[i], R[i], rotation );
				fprintf( out_fd, ">%s\n", seq_id[i] );
				fprintf( out_fd, "%s\n", rotation );
				free ( rotation );
				succ++;
			}
			else
			{
				if ( outliers_filename != NULL )
				{
					if ( unsucc == 0 )
					{
						if ( ! ( outl_fd = fopen ( outliers_filename, "w") ) )
						{
							fprintf ( stderr, " Error: Cannot open file %s!\n", outliers_filename );
							return ( 1 );
						}
					}	
					fprintf( outl_fd, ">%s\n", ( char * ) seq_id[i] );
					fprintf( outl_fd, "%s\n", ( char * ) seq[i] );	
					unsucc++;
				}
			}
				
		}

		fprintf ( stderr, "%d/%d sequences were succesfully rotated.\n", succ, num_seqs );

		if ( fclose ( out_fd ) )
		{
			fprintf( stderr, " Error: file close error!\n");
			return ( 1 );
		}

		if ( outliers_filename != NULL && unsucc > 0 )
		{
			if ( fclose ( outl_fd ) )
			{
				fprintf( stderr, " Error: file close error!\n");
				return ( 1 );
			}
		}

		free ( R );
		for ( i = 0; i < num_seqs; i ++ )
		{
			free ( D[i] );
		}
		free ( D );

	}

	double end = gettime();

        fprintf( stderr, "Elapsed time for processing %d sequence(s): %lf secs.\n", num_seqs, ( end - start ) );
	
	/* Deallocate */
	for ( i = 0; i < num_seqs; i ++ )
	{
		free ( seq[i] );
		free ( seq_id[i] );
	}	
	free ( seq );
	free ( seq_id );
	free ( pat );
        free ( sw . patterns_filename );
	if ( outliers_filename != NULL )
        	free ( sw . outliers_filename );
	if ( text_filename != NULL )
        	free ( sw . text_filename );
	if ( rotations_filename != NULL )
        	free ( sw . rotations_filename );
        free ( sw . output_filename );
        free ( sw . alphabet );

	return ( 0 );
}
