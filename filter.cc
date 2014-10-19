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

#include <math.h>
#include <stdlib.h>
#include "trie.h"						//include header for trie

void fragments ( int i, int j, int f, unsigned int m, int * mf, int * ind )
{
	unsigned int modulo = m % f;
	double nf = ( double ) ( m ) / f;
	int first;
	int last;

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


unsigned int extract_dups ( unsigned char ** xx, unsigned int d, unsigned int * m, unsigned int f, int * mf, int * ind, int * dups )
{
	AlphaMap *      alphabet = NULL;
        Trie *          trie = NULL;

        /* Create an empty alphabet */
        alphabet = alpha_map_new();

        /* Define the alphabet's range */
        alpha_map_add_range ( alphabet, 0, 127 );

        /* Create an empty trie based on the alphabet */
        trie = trie_new ( alphabet );

	unsigned int uniq = 0;

        for( int i = 0; i < d; i++ )
        {
		for ( int j = 0; j < f; j++ )
		{
			TrieData data;
			unsigned int f_id = i * f + j;
			int Nf   = mf[f_id];

			AlphaChar * pf;    
			pf = ( AlphaChar * ) calloc ( ( Nf + 1 ) , sizeof( AlphaChar ) );
			for ( int k = 0; k < Nf; k++ )
			{
				pf[k] = ( AlphaChar ) xx[i][ ind[f_id] + k ];
			}
			pf[Nf] = '\0';

			if ( trie_retrieve ( trie, pf, &data ) != TRUE )
			{
				trie_store ( trie, pf, f_id );
				dups[ f_id ] = -1;		//first occurrence
				uniq ++;
			}
			else
				dups[ f_id ] = data;
	 
			free ( pf );
		}
	}
	
	trie_free ( trie );
        trie = NULL;
        alpha_map_free ( alphabet );
        alphabet = NULL;
	return ( uniq );
}
