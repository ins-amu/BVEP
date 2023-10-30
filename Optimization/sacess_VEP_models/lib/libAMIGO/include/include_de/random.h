//------Constants for the Mersenne Twister---------------------------------
// Period parameters 
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df   // constant vector a
#define UPPER_MASK 0x80000000 // most significant w-r bits
#define LOWER_MASK 0x7fffffff // least significant r bits

// Tempering parameters 
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned long gla_mt[N]; // the array for the state vector  
static int gl_mti=N+1; // gl_mti==N+1 means gla_mt[N] is not initialized 

//-----The Mersenne Twister and its initialization-------------------------
/* A C-program for MT19937: Real number version                */
/*   genrand() generates one pseudorandom real number (double) */
/* which is uniformly distributed on [0,1]-interval, for each  */
/* call. sgenrand(i_seed) set initial values to the working area */
/* of 624 words. Before genrand(), sgenrand(i_seed) must be      */
/* called once. (i_seed is any 32-bit integer except for 0).     */
/* Integer generator is obtained by modifying two lines.       */
/*   Coded by Takuji Nishimura, considering the suggestions by */
/* Topher Cooper and Marc Rieffel in July-Aug. 1997.           */

/* This library is free software; you can redistribute it and/or   */
/* modify it under the terms of the GNU Library General Public     */
/* License as published by the Free Software Foundation; either    */
/* version 2 of the License, or (at your option) any later         */
/* version.                                                        */
/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of  */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            */
/* See the GNU Library General Public License for more details.    */
/* You should have received a copy of the GNU Library General      */
/* Public License along with this library; if not, write to the    */
/* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA   */ 
/* 02111-1307  USA                                                 */

/* Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.       */
/* When you use this, send an email to: matumoto@math.keio.ac.jp   */
/* with an appropriate reference to your work.                     */

// initializing the array with a NONZERO seed 
void sgenrand(unsigned long seed) 
{
    // setting initial seeds to gla_mt[N] using     
    // the generator Line 25 of Table 1 in          
    // [KNUTH 1981, The Art of Computer Programming 
    //    Vol. 2 (2nd Ed.), pp102]                  
    gla_mt[0]= seed & 0xffffffff;
    for (gl_mti=1; gl_mti<N; gl_mti++)
        gla_mt[gl_mti] = (69069 * gla_mt[gl_mti-1]) & 0xffffffff;
}


float genrand()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    // mag01[x] = x * MATRIX_A  for x=0,1 

    if (gl_mti >= N) { // generate N words at one time 
        int kk;

        if (gl_mti == N+1)   // if sgenrand() has not been called,
            sgenrand(4357); // a default initial seed is used   

        for (kk=0;kk<N-M;kk++) {
            y = (gla_mt[kk]&UPPER_MASK)|(gla_mt[kk+1]&LOWER_MASK);
            gla_mt[kk] = gla_mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N-1;kk++) {
            y = (gla_mt[kk]&UPPER_MASK)|(gla_mt[kk+1]&LOWER_MASK);
            gla_mt[kk] = gla_mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (gla_mt[N-1]&UPPER_MASK)|(gla_mt[0]&LOWER_MASK);
        gla_mt[N-1] = gla_mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

        gl_mti = 0;
    }
  
    y = gla_mt[gl_mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

    return ( (float)y / (unsigned long)0xffffffff ); // reals
    /* return y; */ /* for integer generation */
}
