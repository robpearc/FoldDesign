///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "CommonPara.h"

/* setting initial seeds to mt[N] using
 * the generator Line 25 of Table 1 in
 * [KNUTH 1981, The Art of Computer Programming
 * Vol. 2 (2nd Ed.), pp102] */
void sgenrand(
  unsigned long seed
)
{
  mt[0]= seed & 0xffffffff;
  for(mti=1; mti<N; mti++)
    mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}

/* returns a pseudo-random real number uniformly distributed 
 * between 0.0 and 1.0. */
double Random()
{
  const long Q = MODULUS / MULTIPLIER;
  const long R = MODULUS % MULTIPLIER;
  long t;
    
  t = MULTIPLIER * (seed[stream] % Q) - R * (seed[stream] / Q);
  if(t > 0) seed[stream] = t;
  else seed[stream] = t + MODULUS;
  return ((double) seed[stream] / MODULUS);
}

/* Use this function to set the state of all the random number generator 
 * streams by "planting" a sequence of states (seeds), one per stream, 
 * with all states dictated by the state of the default stream. 
 * The sequence of planted states is separated one from the next by 
 * 8,367,782 calls to Random().  */
void PlantSeeds(
  long x
)
{
  const long Q = MODULUS / A256;
  const long R = MODULUS % A256;
  int  j;
  int  s;
    
  initialized = 1;
  s = stream;                   /* remember the current stream */
  SelectStream(0);              /* change to stream 0          */
  PutSeed(x);                   /* set seed[0]                 */
  stream = s;                   /* reset the current stream    */
  for (j = 1; j < STREAMS; j++)
  {
    x = A256 * (seed[j - 1] % Q) - R * (seed[j - 1] / Q);
    if(x > 0) seed[j] = x;
    else seed[j] = x + MODULUS;
  }
}

/* Use this function to set the state of the current random number 
 * generator stream according to the following conventions:
 *    if x > 0 then x is the state (unless too large)
 *    if x < 0 then the state is obtained from the system clock
 *    if x = 0 then the state is to be supplied interactively */
void PutSeed(
  long x
)
{
  char ok = 0;
  if(x > 0) x = x % MODULUS;    /* correct if x is too large  */
  if(x < 0) x = ((unsigned long) time((time_t *) NULL)) % MODULUS;
  if(x == 0)
  {
    while(!ok) 
    {
      printf("\nEnter a positive integer seed (9 digits or less) >> ");
      scanf("%ld", &x);
      ok = (0 < x) && (x < MODULUS);
      if(!ok) printf("\nInput out of range ... try again\n");
    }
  }
  seed[stream] = x;
}

/* Use this function to set the current random number generator
 * stream -- that stream from which the next random number will come. */
void SelectStream(
  int index
)
{
  stream = ((unsigned int) index) % STREAMS;
  if((initialized == 0) && (stream != 0))   /* protect against        */
    PlantSeeds(DEFAULT);                   /* un-initialized streams */
}
