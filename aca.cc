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
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "ahocorasick.h"
#include "beardefs.h"
#include "globals.h"

/****************************************************************************/

//*** 1. Define a call-back function of the MATCH_CALBACK_t type
int match_handler( AC_MATCH_t * m, int * param )
{
	unsigned int j;
	//printf ("@ %ld : ", m->position - 1);

	for ( j = 0; j < m->match_num; j++ )
	{
		//fprintf(stderr, "pos:%d id:%d.\n",  m -> position - m->patterns[j].length, m->patterns[j].rep.number);
		//  ending position of the fragment in text  minus the length of fragment = starting position
		int pos = m -> position - m->patterns[j].length;
		if ( pos > 0 )
		{
			if ( gMatches >= gMax_alloc_matches )
			{
				( gF ) = ( int * ) realloc ( ( gF ),   ( gMax_alloc_matches + ALLOC_SIZE ) * sizeof ( int ) );
				if ( ( gF ) == NULL )
				{
					fprintf ( stderr, " Error: Cannot allocate memory for gF.\n" );
					exit( EXIT_FAILURE );
				}
				( gP ) = ( int * ) realloc ( ( gP ),   ( gMax_alloc_matches + ALLOC_SIZE ) * sizeof ( int ) );
				if ( ( gP ) == NULL )
				{
					fprintf ( stderr, " Error: Cannot allocate memory for gP.\n" );
					exit( EXIT_FAILURE );
				}
				gMax_alloc_matches += ALLOC_SIZE;
			}
			( gF )[ gMatches ] = m->patterns[j].rep.number; 
			( gP )[ gMatches ] = pos; 
			gMatches ++;
		}
	}

	/* to find all matches always return 0 */
	return 0;
	/* return 0 : continue searching
	 * return none zero : stop searching
	 * as soon as you get enough from search results, you can stop search and
	 * return from ac_automata_search() and continue the rest of your program.
	 * e.g. if you only need first N matches, define a counter and return none
	 * zero after the counter exceeds N.
	**/
}


/****************************************************************************/

int filtering ( char * t, unsigned int n, char ** seqs, unsigned int f )
{
	unsigned int i;

	//*** 2. Define AC variables: AC_AUTOMATA_t *, AC_PATTERN_t, and AC_TEXT_t
	AC_AUTOMATA_t * acap;
	AC_PATTERN_t tmp_patt;
	AC_TEXT_t tmp_text;

	//*** 3. Get a new automata
	acap = ac_automata_init ( match_handler );

	for( int i = 0; i < f; ++i )
	{
		if ( strlen( seqs[i] ) > 0 )
		{
			tmp_patt.astring = seqs[i];
			tmp_patt.length = strlen( seqs[i] );
			tmp_patt.rep.number = i; // the id of the fragment
			ac_automata_add (acap, &tmp_patt);
		}
	}

	//*** 5. Finalize automata.
	ac_automata_finalize (acap);

	//*** 6. Set input text
	tmp_text.astring = t;
	tmp_text.length = n;

	//*** 7. Do search
	ac_automata_search (acap, &tmp_text, NULL);
	/* here we pass 0 to our callback function.
	 * if there are variables to pass to call-back function,
	 * you can define a struct that enclose those variables and
	 * send the pointer of the struct as a parameter.
	**/

	//*** 9. Release automata
	ac_automata_release (acap);
	/* if you have finished with the automata and will not use it in the rest
	 * of your program please release it.
	**/

	return 0;
}
