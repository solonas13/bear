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

#ifndef __BEAR__
#define __BEAR__

#define ALLOC_SIZE               10000
#define DEL                     '$'
#define DEL_STR                 "$"
#define ERR                      24
#define PROT                    "ARNDCQEGHILKMFPSTWYVBZX*"   //Proteins alphabet
#define DNA                     "ATGCSWRYKMBVHDN"            //IUPAC alphabet
#define RNA                     "AUGCN"                      //RNA alphabet
#define NUC_SCORING_MATRIX_SIZE 15
#define PRO_SCORING_MATRIX_SIZE 24
#define WORD_LEN 		64
#define max(a,b) ((a) > (b)) ? (a) : (b)
#define min(a,b) ((a) < (b)) ? (a) : (b)
#define nuc_delta(a,b) ((a) == DEL | (b) == DEL) ? 0 : EDNAFULL_matrix[ EDNA[(int)(a)] ][ EDNA[(int)(b)] ]
#define pro_delta(a,b) ((a) == DEL | (b) == DEL) ? 0 : EBLOSUM62_matrix[ BLOSUM[(int)(a)] ][ BLOSUM[(int)(b)] ]

#include <tuple>

typedef unsigned long int WORD;

struct TSwitch
 {
   char               * alphabet;
   char               * patterns_filename;
   char               * text_filename;
   char               * output_filename;
   char               * rotations_filename;
   char               * outliers_filename;
   unsigned int         matrix;
   unsigned int         k;
   unsigned int         D;
   unsigned int         A;
   unsigned int         d;
   double	        O;
   double	        E;
   double	        R;
   double	        P;
   double	        min_sim;
   unsigned int         w;
   unsigned int         T;
   unsigned int         b;
   unsigned int         q;
 };

struct TPat
 {
   unsigned int         start;
   unsigned int		offset;
 };

struct TPOcc
 {
   unsigned int         pos;
   double	        err;
   unsigned int		rot;
 };

typedef std::tuple<int,int,int> mytuple;

typedef int32_t INT;

unsigned int LCParray ( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP );
void partitioning ( INT i, INT j, INT f, INT m, INT * mf, INT * ind );
unsigned int circular_sequence_comparison (  unsigned char * x, unsigned char * y, struct TSwitch  sw, unsigned int * rotation, unsigned int * distance );
unsigned int sacsc_refinement (  unsigned char * x, unsigned char * y, struct TSwitch  sw, unsigned int * rotation, unsigned int * distance );
unsigned int nw_refine ( unsigned char * p, unsigned int m, unsigned char * t, unsigned int n, double o, double e, double * score, struct TSwitch sw );
int refine ( unsigned char * x, unsigned int m, unsigned char * y, unsigned int n, struct TSwitch sw );

unsigned int pcsa_ed( unsigned char * x, unsigned char * y, struct TSwitch  sw, unsigned int * rotation, unsigned int * distance );
unsigned int pcsa_hd( unsigned char * x, unsigned char * y, struct TSwitch  sw, unsigned int * rotation, unsigned int * distance );
unsigned int macsmf_ed( unsigned char ** x, unsigned char * t, struct TSwitch  sw, TPOcc *** POcc, unsigned int ** NOcc );
unsigned int macsmf_hd( unsigned char ** x, unsigned char * t, struct TSwitch  sw, TPOcc *** POcc, unsigned int ** NOcc );
unsigned int cyc_nw_ls ( unsigned char * x, unsigned int m, unsigned char * y, unsigned int n, struct TSwitch sw, double * score, int * rot );
unsigned int sw_ls ( unsigned char * p, unsigned int m, unsigned  char * t, unsigned int n, struct TSwitch sw, unsigned int * ii, unsigned int * jj, double * similarity );
unsigned int hamming_distance ( unsigned char * x, unsigned int m, unsigned char * y,  unsigned int n, struct TSwitch  sw, unsigned int * S );
unsigned int edit_distance ( unsigned char * x, unsigned int m, unsigned char * y, unsigned int n, struct TSwitch  sw, unsigned int ** D );
unsigned int edit_distance_wbt ( unsigned char * x, unsigned int m, unsigned char * y, unsigned int n, struct TSwitch  sw, unsigned int ** D, int ** H );
unsigned int r_errors_vec	( unsigned int ** D, unsigned int m, unsigned int n, struct TSwitch  sw, int * Sl );
unsigned int l_errors_vec	( unsigned int ** D, unsigned int m, unsigned int n, int ** H, struct TSwitch  sw, int * Sr, int * Pr );
unsigned int bcf_maxshift_hd_ls ( unsigned char * p, unsigned int m, unsigned  char * t, unsigned int n, unsigned int l, unsigned int * ii, unsigned int * jj, unsigned int * distance );
unsigned int bcf_maxshift_ed_ls ( unsigned char * p, unsigned int m, unsigned  char * t, unsigned int n, unsigned int l, unsigned int * ii, unsigned int * jj, unsigned int * distance );
unsigned int upgma_dist ( TPOcc ** POcc, unsigned int d, struct TSwitch  sw, int * R, unsigned char ** seq );
unsigned int upgma_sim ( TPOcc ** POcc, unsigned int d, struct TSwitch  sw, int * R, unsigned char ** seq );
//int nuc_delta ( char a, char b );
//int pro_delta ( char a, char b );
void init_substitution_score_tables ( void );
void create_rotation ( unsigned char * x, unsigned int offset, unsigned char * rotation );
void create_backward_rotation ( unsigned char * x, unsigned int offset, unsigned char * rotation );
TPOcc * unique ( TPOcc * first, TPOcc * last );
int decode_switches ( int argc, char * argv [], struct TSwitch * sw );
double gettime ( void );
void usage ( void );

#endif

