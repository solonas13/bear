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
#include "EDNAFULL.h"
#include "EBLOSUM62.h"

unsigned int EDNA[90];
unsigned int BLOSUM[91];

void init_substitution_score_tables ( void )
{
    unsigned int i;
    char edna[] = DNA;
    for ( i = 0; i < NUC_SCORING_MATRIX_SIZE; i ++ ) {
	EDNA[(int)edna[i]] = i;
    }
    EDNA[(int)'U'] = 1; //Setting RNA U = T
    char blosum[] = PROT;
    for ( i = 0; i < PRO_SCORING_MATRIX_SIZE; i ++ ) {
	BLOSUM[(int)blosum[i]] = i;
    }
}

/*
 * The new methods have any base that matches a DEL character return an exact
 * matching score for the first nucleotide base or amino acid (A).
 */
/* Returns the score for matching character a and b based on EDNAFULL matrix */
int nuc_delta ( char a, char b )
{
    if ( a == DEL || b == DEL ) {
	a = DNA[0];
	b = DNA[0];
    }
    return EDNAFULL_matrix[ EDNA[(int)a] ][ EDNA[(int)b] ];
}

/* Returns the score for matching character a and b based on EBLOSUM62 matrix */
int pro_delta ( char a, char b )
{
    if ( a == DEL || b == DEL ) {
	a = PROT[0];
	b = PROT[0];
    }
    return EBLOSUM62_matrix[ BLOSUM[(int)a] ][ BLOSUM[(int)b] ];
}

///*
// * This is the old code where an exact match of DEL chars returned a very high
// * score and a mismatch with DEL a very bad score, and was used to enforce
// * excellent matches in the middle of strings and small mismatches at either
// * end of the strings; used in the refine() function
// */
// /* Returns the score for matching character a and b based on EDNAFULL matrix */
// int nuc_delta ( char a, char b )
//  {
//     if ( a == DEL && b == DEL ) {
// 	return 10;
//     }
//     else if ( a == DEL || b == DEL ) {
// 	return -10;
//     }
//     return EDNAFULL_matrix[ EDNA[(int)a] ][ EDNA[(int)b] ];
//  }
// 
//  /* Returns the score for matching character a and b based on EBLOSUM62 matrix */
// int pro_delta ( char a, char b )
// {
//     if ( a == DEL && b == DEL ) {
// 	return 10;
//     }
//     else if ( a == DEL || b == DEL ) {
// 	return -10;
//     }
//     return EBLOSUM62_matrix[ BLOSUM[(int)a] ][ BLOSUM[(int)b] ];
// }
