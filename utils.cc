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
#include <sys/time.h>
#include "beardefs.h"

double gettime( void )
{
    struct timeval ttime;
    gettimeofday( &ttime , 0 );
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
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

void create_rotation ( unsigned char * x, unsigned int offset, unsigned char * rotation )
{
	unsigned int m = strlen ( ( char * ) x );
	memmove ( &rotation[0], &x[offset], m - offset );
	memmove ( &rotation[m - offset], &x[0], offset );
	rotation[m] = '\0';
}

void create_backward_rotation ( unsigned char * x, unsigned int offset, unsigned char * rotation )
{
	unsigned int m = strlen ( ( char * ) x );
	memmove ( &rotation[0], &x[m - offset], offset );
	memmove ( &rotation[offset], &x[0], m - offset );
	rotation[m] = '\0';
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

