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

#define ALLOC_SIZE               10000
#define DEL                     '$'
#define DEL_STR                 "$"
#define ERR                      24
#define PROT                    "ARNDCQEGHILKMFPSTWYVBZX*U"   //Proteins alphabet
#define DNA                     "ATGCSWRYKMBVHDN"            //IUPAC alphabet
#define NUC_SCORING_MATRIX_SIZE 15
#define PRO_SCORING_MATRIX_SIZE 24
#define WORD_LEN 		64
#define max(a,b) ((a) > (b)) ? (a) : (b)
#define min(a,b) ((a) < (b)) ? (a) : (b)

#include <tuple>

using namespace std;
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
   double         	O;
   double	        E;
   double         	R;
   double         	min_sim;
   unsigned int         w;
   unsigned int		T;	
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

typedef tuple<int,int,int> mytuple;

unsigned int macsmf_ed( unsigned char ** x, unsigned char * t, struct TSwitch  sw, TPOcc *** POcc, unsigned int ** NOcc );
unsigned int macsmf_hd( unsigned char ** x, unsigned char * t, struct TSwitch  sw, TPOcc *** POcc, unsigned int ** NOcc );
unsigned int nw ( unsigned char * x, unsigned int m, unsigned char * y, unsigned int n, struct TSwitch  sw, double ** D );
unsigned int nw_wbt ( unsigned char * x, unsigned int m, unsigned char * y, unsigned int n, struct TSwitch  sw, double ** D, int ** H);
unsigned int cyc_nw_ls ( unsigned char * x, unsigned int m, unsigned char * y, unsigned int n, struct TSwitch sw, double * score, int * rot );
unsigned int sw ( unsigned char * p, unsigned int m, unsigned  char * t, unsigned int n, struct TSwitch sw, TPOcc * D );
unsigned int sw_ls ( unsigned char * p, unsigned int m, unsigned  char * t, unsigned int n, struct TSwitch sw, TPOcc * D );
unsigned int hamming_distance ( unsigned char * x, unsigned int m, unsigned char * y,  unsigned int n, struct TSwitch  sw, unsigned int * S );
unsigned int edit_distance ( unsigned char * x, unsigned int m, unsigned char * y, unsigned int n, struct TSwitch  sw, unsigned int ** D );
unsigned int edit_distance_wbt ( unsigned char * x, unsigned int m, unsigned char * y, unsigned int n, struct TSwitch  sw, unsigned int ** D, int ** H );
unsigned int r_errors_vec	( unsigned int ** D, unsigned int m, unsigned int n, struct TSwitch  sw, int * Sl );
unsigned int l_errors_vec	( unsigned int ** D, unsigned int m, unsigned int n, int ** H, struct TSwitch  sw, int * Sr, int * Pr );
unsigned int upgma_dist2 ( TPOcc ** POcc, unsigned int d, struct TSwitch  sw, int * R );
unsigned int upgma_dist ( TPOcc ** POcc, unsigned int d, struct TSwitch  sw, int * R, unsigned char ** seq );
unsigned int upgma_sim2 ( TPOcc ** POcc, unsigned int d, struct TSwitch  sw, int * R );
unsigned int upgma_sim ( TPOcc ** POcc, unsigned int d, struct TSwitch  sw, int * R, unsigned char ** seq );
unsigned int maxshift_hd ( unsigned char * p, unsigned int m, unsigned  char * t, unsigned int n, struct TSwitch sw, TPOcc * D );
unsigned int maxshift_ed ( unsigned char * p, unsigned int m, unsigned  char * t, unsigned int n, struct TSwitch sw, TPOcc * D );
unsigned int bcf_maxshift_ed ( unsigned char * p, unsigned int m, unsigned  char * t, unsigned int n, struct TSwitch sw, TPOcc * D );
unsigned int bcf_maxshift_hd ( unsigned char * p, unsigned int m, unsigned  char * t, unsigned int n, struct TSwitch sw, TPOcc * D );
unsigned int bcf_maxshift_hd_ls ( unsigned char * p, unsigned int m, unsigned  char * t, unsigned int n, struct TSwitch sw, TPOcc * D );
unsigned int bcf_maxshift_ed_ls ( unsigned char * p, unsigned int m, unsigned  char * t, unsigned int n, struct TSwitch sw, TPOcc * D );
TPOcc * unique ( TPOcc * first, TPOcc * last );
int binSearch( int value, TPat * array, int num);
unsigned int nuc_char_to_index ( char a );
unsigned int pro_char_to_index ( char a );
int nuc_delta ( char a, char b );
int pro_delta ( char a, char b );
unsigned int create_rotation ( unsigned char * x, unsigned int offset, unsigned char * rotation );
int decode_switches ( int argc, char * argv [], struct TSwitch * sw );
double gettime( void );
void usage ( void );
