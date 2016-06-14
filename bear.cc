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
#include <libflasm.h>
#include "globals.h"
#include "EDNAFULL.h"
#include "EBLOSUM62.h"
#include "beardefs.h"

using namespace std;
using namespace libflasm;

int main(int argc, char **argv)
{
	struct TSwitch  sw;

	FILE *          in_fd;                  // the input file descriptor
	FILE *          out_fd;                 // the output file descriptor
	FILE *          outl_fd;                // the outliers file descriptor
	//FILE *          rotsAndScores_fd;       // AR used to store info for debugging
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
	unsigned int    d, q, b;	
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
                else if ( ! strcmp ( "RNA", sw . alphabet ) )   { alphabet = ( char * ) RNA;  sw . matrix = 0; }
                else if ( ! strcmp ( "PROT", sw . alphabet ) )  { alphabet = ( char * ) PROT; sw . matrix = 1; }
                else
                {
                        fprintf ( stderr, " Error: alphabet argument a should be `DNA' for nucleotide sequences or `PROT' for protein sequences or `USR' for sequences over a user-defined alphabet!\n" );
                        return ( 1 );
                }

		d       = sw . d;
		if ( d > 3 )	
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

		if ( d == 3 )
		{
			if ( sw . P < 0 || sw . P >= 49.0 )
			{
				fprintf ( stderr, " Error: The optional percentage flag should be in the range of 0 to 49%%.\n" );
				return ( 1 );
			}
			if ( sw . q < 2 || sw . q >= sw . b )
			{
				fprintf ( stderr, " Error: The length of the q-gram must be reasonable.\n" );
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
		if ( sw . D == 0 ) //AR Edit to if ( d==0 )
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
		/* D[i,j] stores distance( x^{rot}_{i}, x_{j} ) = err */
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
						unsigned int distance = ( int ) DBL_MAX;
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
						unsigned int distance = ( int ) DBL_MAX;
						unsigned int rotation = 0;

						pcsa_ed ( seq[i], seq[j], sw, &rotation, &distance );
						
						D[i][j] . err = distance;
						D[i][j] . rot = rotation;		
					}
				}
			}
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

				unsigned char * xiw;
				xiw = ( unsigned char * ) calloc (  m + sw . w, sizeof( unsigned char ) );
				memmove( &xiw[0], seq[i], m );
				memmove( &xiw[m], seq[i], sw . w - 1 );
				xiw[m + sw . w - 1] = '\0';
				unsigned int mm = m + sw . w - 1;

				for ( int j = 0; j < num_seqs; j ++ )
				{
					if ( i == j ) continue;

					/* Here we compute these coordinates using libFLASM and store it to M[ii] */
					unsigned int n = strlen ( ( char * ) seq[j] );

					D[i][j] . err = m - 1;
					unsigned int ii, jj;
					unsigned int distance = ( unsigned int ) DBL_MAX;
					unsigned int l = min ( sw . w, m );

					if ( sw . D == 0 )
					{
						ResultTupleSet results = flasm_hd ( xiw, mm, seq[j], n, l, sw . k, false );
						if ( results . size () > 0 ) {
							ResultTupleSetIterator it = results . begin();
							ii = ( * it ) . pos_t;
							jj = ( * it ) . pos_x;
							distance = ( * it ) . error;
						}
					}
					else
					{
						ResultTupleSet results = flasm_ed ( xiw, mm, seq[j], n, l, sw . k, false );
						if ( results . size () > 0 ) {
							ResultTupleSetIterator it = results . begin();
							ii = ( * it ) . pos_t;
							jj = ( * it ) . pos_x;
							distance = ( * it ) . error;
						}
					}

					D[i][j] . err = distance;

					/* If ii >= jj then it is easy to find the rotation: ii - jj */
					if ( ii >= jj )
					{
						D[i][j] . rot = ii - jj;
					}
					/* Otherwise we would need to shift jj - ii or as many letters as we have on the right of position i - sw. w + 1 */
					else
					{
						int a = jj - ii;
						int b = m - ii - 1;
						int c = min ( a, b );
						D[i][j] . rot = m - c;
					}
				}
				free ( xiw );
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

					unsigned char * xx;
					xx = ( unsigned char * ) calloc (  m + m, sizeof( unsigned char ) );
					memmove( &xx[0], seq[i], m );
					memmove( &xx[m], seq[i], m - 1 );
					xx[m + m - 1] = '\0';
					unsigned int mm = m + m - 1;

					for ( int j = 0; j < num_seqs; j ++ )
					{
						if ( i == j ) continue;
						unsigned int n = strlen ( ( char * ) seq[j] );

						D[i][j] . err = -DBL_MAX;
						unsigned int ii, jj;
						double similarity = -DBL_MAX;

						sw_ls ( xx, mm, seq[j], n, sw, &ii, &jj, &similarity );

						D[i][j] . err = similarity;

						/* If ii >= jj then it is easy to find the rotation: ii - jj */
						if ( ii >= jj )
						{
							D[i][j] . rot = ii - jj;
						}
						/* Otherwise we would need to shift jj - ii or as many letters as we have on the right of position i - sw. w + 1 */
						else
						{
							int a = jj - ii;
							int b = m - ii - 1;
							int c = min ( a, b );
							D[i][j] . rot = m - c;		
						}

					}
					free ( xx );
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
				
				//AR open file to write debugging info
				/*if ( ! ( rotsAndScores_fd = fopen ( "rotsAndScores.2.txt", "w") ) )
				{
					fprintf ( stderr, " Error: Cannot open file rotsAndScores.2.txt!\n" );
					return ( 1 );
				}*/

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
						D[i][j] . rot = rot;
						
						//AR write out scores and rotations
						//fprintf ( rotsAndScores_fd, "i: %02d, j: %02d, dst: %08.2f, rot: %04u\n", i, j, D[i][j] . err, D[i][j] . rot );
					}
				}

				//AR close debugging file
				//fclose(rotsAndScores_fd);
			}	
		}

		if ( d == 3 )
		{
			fprintf ( stderr, " Comparison model: blockwise q-gram distance.\n" );
			if ( ( D = ( TPOcc ** ) calloc ( ( num_seqs ) , sizeof( TPOcc * ) ) ) == NULL )
			{
				fprintf( stderr, " Error: Cannot allocate memory!\n" );
				exit( EXIT_FAILURE );
			}

			fprintf ( stderr, " Length of block = %d, q = %d and refinement %1.2f%%.\n", sw . b, sw . q, sw . P );

			//AR open file to write debugging info
			/*if ( ! ( rotsAndScores_fd = fopen ( "rotsAndScores.3.txt", "w") ) )
			{
				fprintf ( stderr, " Error: Cannot open file rotsAndScores.3.txt!\n" );
				return ( 1 );
			}*/
			
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

					unsigned int n = strlen ( ( char * ) seq[j] );
					if ( sw . b > m - sw . q + 1  || sw . b > n - sw . q + 1 )
					{
						fprintf( stderr, " Error: Illegal block length.\n" );
						exit ( 1 );
					}

					sw . k = n + m;
					unsigned int distance = ( int ) DBL_MAX;
					unsigned int rotation = 0;

					circular_sequence_comparison ( seq[i], seq[j], sw, &rotation, &distance );

					D[i][j] . err = (double) distance;
					D[i][j] . rot = rotation;

					//AR write out scores and rotations
					//fprintf ( rotsAndScores_fd, "i: %02d, j: %02d, dst: %08.2f, rot: %04u\n", i, j, D[i][j] . err, D[i][j] . rot );
				}
			}
			
			//AR close debugging file
			//fclose(rotsAndScores_fd);
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
		if ( d == 0 || d == 1 || d == 3 )
			upgma_dist ( D, num_seqs, sw, R, seq );
		else if ( d == 2 )
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
