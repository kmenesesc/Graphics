/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.
   Converstion to C++ by Sam Buss, 2007.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).
   In Buss's implemention, use the corresponding constructors instead

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   This software is distributed "as is" and carries no warranties.

   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)

   Sam Buss (C++ conversion only)
   sbuss@math.ucsd.edu
*/

#ifndef MERSENNE_TWISTER_H
#define MERSENNE_TWISTER_H

class MersenneTwister
{
public:
	MersenneTwister();
	MersenneTwister( unsigned long seed );
	MersenneTwister( unsigned const long init_key[], int key_length );

	unsigned long Rand32();
	long Rand31();

	// Call these to re-seed the random numbers
	void init_genrand(unsigned long s);
	void init_by_array(unsigned const long init_key[], int key_length);

private:
	/* Period parameters */  
	static const int N;
	static const int M;
	static const unsigned long MATRIX_A;						/* constant vector a */
	static const unsigned long UPPER_MASK;						/* most significant w-r bits */
	static const unsigned long LOWER_MASK;						/* least significant r bits */
	static unsigned long mag01[2]; // = {0x0UL, MATRIX_A};		/* mag01[x] = x * MATRIX_A  for x=0,1 */

private:
	unsigned long mt[624];								/* the array for the state vector : N=624 */
	int mti;											/* mti==N+1 means mt[N] is not initialized */
};

inline MersenneTwister::MersenneTwister()
{
	mti = N+1;												/* mti==N+1 means mt[N] is not initialized */
	init_genrand( 5489UL );									/* a default initial seed is used */
}

inline MersenneTwister::MersenneTwister( unsigned long seed )
{
	mti = N+1;												/* mti==N+1 means mt[N] is not initialized */
	init_genrand( seed );
}

inline MersenneTwister::MersenneTwister( unsigned const long init_key[], int key_length )
{
	mti = N+1;												/* mti==N+1 means mt[N] is not initialized */
	init_by_array( init_key, key_length);
}

/* generates a random number on [0,0xffffffff]-interval */
inline unsigned long MersenneTwister::Rand32()
{
    unsigned long y;

    if (mti >= N) { /* generate N words at one time */
        int kk;
        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
inline long MersenneTwister::Rand31()
{
    return (long)(Rand32()>>1);
}


#endif // MERSENNE_TWISTER_H