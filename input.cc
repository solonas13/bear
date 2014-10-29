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
#include "beardefs.h"

static struct option long_options[] =
 {
   { "alphabet",                required_argument, NULL, 'a' },
   { "seqs-file",               required_argument, NULL, 'p' },
   { "ref-file",                required_argument, NULL, 't' },
   { "output-file",             required_argument, NULL, 'o' },
   { "outliers-file",           required_argument, NULL, 'l' },
   { "max-dist",                required_argument, NULL, 'k' },
   { "max-gap",                 required_argument, NULL, 'g' },
   { "fac-len",                 required_argument, NULL, 'w' },
   { "sim-rat",                 required_argument, NULL, 'R' },
   { "opn-gap",                 required_argument, NULL, 'O' },
   { "ext-gap",                 required_argument, NULL, 'E' },
   { "rot-file",                required_argument, NULL, 'r' },
   { "num-threads",             required_argument, NULL, 'T' },
   { "help",                    no_argument,       NULL, 'h' },
   { NULL,                      0,                 NULL, 0   }
 };


/* 
Decode the input switches 
*/
int decode_switches ( int argc, char * argv [], struct TSwitch * sw )
 {
   int          oi;
   int          opt;
   double       val;
   char       * ep;
   int          args;

   /* initialisation */
   sw -> alphabet                       = NULL;
   sw -> patterns_filename              = NULL;
   sw -> text_filename                  = NULL;
   sw -> output_filename                = NULL;
   sw -> rotations_filename             = NULL;
   sw -> outliers_filename              = NULL;
   sw -> gaps                           = 2;
   sw -> T                              = 1;
   sw -> k                              = 0;
   sw -> w                              = 45;
   sw -> O                              = -10.0;
   sw -> E                              = -1.0;
   sw -> R                              = 0.25;
   sw -> d                              = 0;
   args = 0;

   while ( ( opt = getopt_long ( argc, argv, "a:p:t:o:k:T:g:w:r:d:l:O:E:R:h", long_options, &oi ) ) != - 1 )
    {
      switch ( opt )
       {
         case 'a':
           sw -> alphabet = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> alphabet, optarg );
           args ++;
           break;

         case 'p':
           sw -> patterns_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> patterns_filename, optarg );
           args ++;
           break;

         case 'o':
           sw -> output_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> output_filename, optarg );
           args ++;
          break;

          case 'k':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> k = val;
           args ++;
           break;

         case 'd':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> d = val;
           args ++;
           break;

         case 'g':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> gaps = val;
           break;

         case 'w':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> w = val;
           break;

         case 'O':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> O = val;
           break;

         case 'E':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> E = val;
           break;

         case 'R':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> R = val;
           break;

          case 'T':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> T = val;
           break;

         case 't':
           sw -> text_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> text_filename, optarg );
           break;

         case 'r':
           sw -> rotations_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> rotations_filename, optarg );
           break;

         case 'l':
           sw -> outliers_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> outliers_filename, optarg );
           break;

         case 'h':
           return ( 0 );
       }
    }

   if ( args < 5 )
     {
       usage ();
       exit ( 1 );
     }
   else
     return ( optind );
 }


/* 
Usage of the tool 
*/
void usage ( void )
 {
   fprintf ( stdout, " Usage: bear <options>\n" );
   fprintf ( stdout, " Standard (Mandatory):\n" );
   fprintf ( stdout, "  -a, --alphabet            <str>     `DNA' for nucleotide  sequences or `PROT' for\n"
                     "                                      protein  sequences. \n" );
   fprintf ( stdout, "  -p, --seqs-file           <str>     MultiFASTA input filename.\n" );
   fprintf ( stdout, "  -o, --output-file         <str>     Output filename.\n" );
   fprintf ( stdout, "  -k, --max-dist            <int>     Maximum distance between pairs of sequences.\n");
   fprintf ( stdout, "  -d, --cmp-mod             <int>     Comparison model (0 for Hamming distance; or\n"
                     "                                      1 for affine gap penalty with sub matrices; or\n"
                     "                                      2 for More Divergent Sequences).\n" );
   fprintf ( stdout, " Optional:\n" );
   fprintf ( stdout, "  -w, --fac-len             <int>     Maximum length of factor to be searched with\n"
                     "                                      option `-d 2' (default: 45). \n" );
   fprintf ( stdout, "  -T, --num-threads         <int>     Number of threads to be used (default: 1).\n" );
   fprintf ( stdout, "  -g, --max-gap             <int>     Maximum total length of gaps between pairs\n"
                     "                                      of sequences (default: 2). \n" );
   fprintf ( stdout, "  -O, --opn-gap             <dbl>     Affine gap open penalty (default: -10.0).\n" );
   fprintf ( stdout, "  -E, --ext-gap             <dbl>     Affine gap extension penalty (default: -1.0).\n" );
   fprintf ( stdout, "  -R, --sim-rat             <dbl>     Ratio of minimum allowed to maximal score\n"
                     "                                      based on the sub matrix used (default: 0.25).\n" );
   fprintf ( stdout, "  -l, --outliers-file       <str>     Outliers filename --- outputs a file with\n"
                     "                                      the sequences which were not rotated. \n" );
   fprintf ( stdout, " Miscellaneous:\n" );
   fprintf ( stdout, "  -t, --ref-file            <str>     Reference sequence filename in FASTA format.\n" );
   fprintf ( stdout, "  -r, --rot-file            <str>     MultiFASTA filename --- outputs a file with\n"
                     "                                      random rotations of the input sequences. \n" );
 }
