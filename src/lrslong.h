/* lrslong.h      (lrs long integer arithmetic library              */
/* Copyright: David Avis 2000, avis@cs.mcgill.ca                    */
/* Version 4.0, February 17, 2000                                   */
/* Version 4.1, November 15, 2020 for bigger  mulint range          */

/* This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Suite 500, Boston, MA  02110-1335, USA.
 */
/******************************************************************************/
/*  See http://cgm.cs.mcgill.ca/~avis/C/lrs.html for lrs usage instructions   */
/******************************************************************************/
/* This package contains the extended precision routines used by lrs
   and some other miscellaneous routines. The maximum precision depends on
   the parameter MAX_DIGITS defined below, with usual default of 255L. This
   gives a maximum of 1020 decimal digits on 32 bit machines. The procedure
   lrs_mp_init(dec_digits) may set a smaller number of dec_digits, and this
   is useful if arrays or matrices will be used.
 */

/***********/
/* defines */
/***********/
/*
   this is number of longwords. Increasing this won't cost you that much
   since only variables other than the A matrix are allocated this size.
   Changing affects running time in small but not very predictable ways.
 */

#define MAX_DIGITS 255L

/*
   this is in decimal digits, you pay in memory if you increase this,
   unless you override by a line with
   digits n
   before the begin line of your file.
 */
#define DEFAULT_DIGITS 100L

/**********MACHINE DEPENDENT CONSTANTS***********/
/* MAXD is 2^(k-1)-1 where k is word size       */
/* MAXDm is sqrt(2^(k-1)-1) where k is word size*/
/* MAXDl is 2^(k/2-1)-1 where k is word size    */
/* MAXDa is 2^(k-2)-1 where k is word size      */
/* MAXD must be at least 2*BASE^2               */
/* If BASE is 10^k, use "%k.ku" for FORMAT      */
/* INTSIZE is number of bytes for integer       */
/* 64/128 bit arithmetic                        */
/************************************************/
#ifdef B128
/* 128 bit machines */ /* compiler does not accept big constants! */
extern __int128 MAXDm, MAXDl,
    MAXDa; /* set correctly in lrs_mp_init in lrslong.c */

/* max power of 10 fitting in signed int64 */
#define P10_INT64 1000000000000000000ULL

#define MAXD 9223372036854775807L /* should be 2^127 -1 but is  2^63 - 1 */
#define BASE 1000000000L
#define FORMAT "%9.9u"
#define BASE_DIG 9
#define INTSIZE 16L
#define BIT "128bit"
#else
/* 64 bit machines */
#define MAXD 9223372036854775807LL  /* 2^63 - 1 */
#define MAXDl 2147483647LL          /* 2^31 - 1 */
#define MAXDm 3037000499LL          /* sqrt(2^63 - 1) */
#define MAXDa 4611686018427387903LL /* 2^62 - 1 */
#define BASE 1000000000L
#define FORMAT "%9.9u"
#define BASE_DIG 9
#define INTSIZE 16L
#define BIT "64bit"
#endif

#define MAXINPUT 1000 /*max length of any input rational */

#define POS 1L
#define NEG -1L
#ifndef TRUE
#define TRUE 1L
#endif
#ifndef FALSE
#define FALSE 0L
#endif
#define ONE 1L
#define TWO 2L
#define ZERO 0L

/**********************************/
/*         MACROS                 */
/* dependent on mp implementation */
/**********************************/

#ifdef SAFE
/* lazy but fast overflow checking */

#define mpsafem(a, b)                                                          \
  *(a) > MAXDm || *(b) > MAXDm || *(a) < -MAXDm || *(b) < -MAXDm
#define mpsafel(a, b)                                                          \
  *(a) > MAXDl || *(b) > MAXDl || *(a) < -MAXDl || *(b) < -MAXDl
#define mpsafea(a, b)                                                          \
  *(a) > MAXDa || *(b) > MAXDa || *(a) < -MAXDa || *(b) < -MAXDa

#ifdef DEBUG
#define mperrorm(a, b)                                                         \
  fprintf(stdout, "  : max(|a|,|b|) > %ld\n", MAXDm);                          \
  lrs_overflow(1)
#define mperrorl(a, b)                                                         \
  fprintf(stdout, "  : max(|a|,|b|) > %ld\n", MAXDl);                          \
  lrs_overflow(1)
#define mperrora(a, b)                                                         \
  fprintf(stdout, "  : max(|a|,|b|) > %ld\n", MAXDa);                          \
  lrs_overflow(1)
#define linint(a, ka, b, kb)                                                   \
  {                                                                            \
    if (mpsafel(a, b)) {                                                       \
      fprintf(stdout, "\n*linint ");                                           \
      mperrorl(a, b);                                                          \
    } else                                                                     \
      *(a) = *(a) * ka + *(b) * kb;                                            \
  }
#define mulint(a, b, c)                                                        \
  {                                                                            \
    if (mpsafem(a, b)) {                                                       \
      fprintf(stdout, "\n*mulint ");                                           \
      mperrorm(a, b);                                                          \
    } else                                                                     \
      *(c) = *(a) * *(b);                                                      \
  }
#define addint(a, b, c)                                                        \
  {                                                                            \
    if (mpsafea(a, b)) {                                                       \
      fprintf(stdout, "\n*addint ");                                           \
      mperrora(a, b);                                                          \
    } else                                                                     \
      *(c) = *(a) + *(b);                                                      \
  }
#define subint(a, b, c)                                                        \
  {                                                                            \
    if (mpsafea(a, b)) {                                                       \
      fprintf(stdout, "\n*subint ");                                           \
      mperrora(a, b);                                                          \
    } else                                                                     \
      *(c) = *(a) - *(b);                                                      \
  }
#define decint(a, b)                                                           \
  {                                                                            \
    if (mpsafea(a, b)) {                                                       \
      fprintf(stdout, "\n*decint ");                                           \
      mperrora(a, b);                                                          \
    } else                                                                     \
      *(a) = *(a) - *(b);                                                      \
  }
#define qpiv(a, b, c, d, e)                                                    \
  {                                                                            \
    if (mpsafel(a, b) || mpsafel(c, d)) {                                      \
      fprintf(stdout, "\n*qpiv ");                                             \
      mperrorl(a, b);                                                          \
      mperrorl(c, d);                                                          \
    }; else                                                                    \
      *(a) = (*(a) * *(b) - *(c) * *(d)) / (*e);                               \
  }
#else
#define qpiv(a, b, c, d, e)                                                    \
  {                                                                            \
    if (mpsafel(a, b) || mpsafel(c, d))                                        \
      lrs_overflow(1);                                                         \
    else                                                                       \
      *(a) = (*(a) * *(b) - *(c) * *(d)) / (*e);                               \
  }
#define linint(a, ka, b, kb)                                                   \
  {                                                                            \
    if (mpsafel(a, b))                                                         \
      lrs_overflow(1);                                                         \
    else                                                                       \
      *(a) = *(a) * ka + *(b) * kb;                                            \
  }
#define mulint(a, b, c)                                                        \
  {                                                                            \
    if (mpsafem(a, b))                                                         \
      lrs_overflow(1);                                                         \
    else                                                                       \
      *(c) = *(a) * *(b);                                                      \
  }
#define addint(a, b, c)                                                        \
  {                                                                            \
    if (mpsafea(a, b))                                                         \
      lrs_overflow(1);                                                         \
    else                                                                       \
      *(c) = *(a) + *(b);                                                      \
  }
#define decint(a, b)                                                           \
  {                                                                            \
    if (mpsafea(a, b))                                                         \
      lrs_overflow(1);                                                         \
    else                                                                       \
      *(a) = *(a) - *(b);                                                      \
  }
#endif

#else
/* unprotected routines */
#define qpiv(a, b, c, d, e) *(a) = (*(a) * *(b) - *(c) * *(d)) / (*e)
#define addint(a, b, c) *(c) = *(a) + *(b)
#define linint(a, ka, b, kb) *(a) = *(a) * ka + *(b) * kb
#define mulint(a, b, c) *(c) = *(a) * *(b)
#endif

// #define unchecked_decint(a, b)  *(a) = *(a) - *(b)  /* v.7.1 only safe if a,b
// come from mulint */
#define unchecked_decint(a, b)                                                 \
  decint(a, b) /* v 7.2 has larger range for mulint       */
#define divint(a, b, c)                                                        \
  {                                                                            \
    *(c) = *(a) / *(b);                                                        \
    *(a) = *(a) % *(b);                                                        \
  }
#define exactdivint(a, b, c) *(c) = *(a) / *(b);

#define abs128(a) (a > 0 ? a : -1 * a)
#define changesign(a) (*(a) = -*(a))
#define copy(a, b) ((a)[0] = (b)[0])
#define mp_greater(a, b) (*(a) > *(b))
#define itomp(in, a) *(a) = in

#define one(a) (*(a) == 1)
#define negative(a) (*(a) < 0)
#define normalize(a) (void)0
#define positive(a) (*(a) > 0)
#define sign(a) (*(a) < 0 ? NEG : POS)
#ifndef B128
#define storesign(a, sa) (*(a) = labs(*(a)) * sa)
#else
#define storesign(a, sa) (*(a) = abs128(*(a)) * sa)
#endif
#define zero(a) (*(a) == 0)

/*
 *  convert between decimal and machine (longword digits). Notice lovely
 *  implementation of ceiling function :-)
 */
#define DEC2DIG(d) ((d) % BASE_DIG ? (d) / BASE_DIG + 1 : (d) / BASE_DIG)
#define DIG2DEC(d) ((d) * BASE_DIG)

#include <stdlib.h> /* labs */

#define CALLOC(n, s) xcalloc(n, s, __LINE__, __FILE__)

/*************/
/* typedefs  */
/*************/
#ifndef B128
typedef long long lrs_mp[1]; /* type lrs_mp holds one long integer */
typedef long long *lrs_mp_t;
typedef long long **lrs_mp_vector;
typedef long long ***lrs_mp_matrix;
#else
typedef __int128 lrs_mp[1]; /* type lrs_mp holds one 128-bit integer */
typedef __int128 *lrs_mp_t;
typedef __int128 **lrs_mp_vector;
typedef __int128 ***lrs_mp_matrix;
#endif

/*********************/
/*global variables   */
/*********************/

extern __thread long lrs_digits;        /* max permitted no. of digits   */
extern __thread long lrs_record_digits; /* this is the biggest acheived so far.     */

/*********************************************************/
/* Initialization and allocation procedures - must use!  */
/******************************************************* */

// void mulint(lrs_mp a, lrs_mp b, lrs_mp c);

long lrs_mp_init(long dec_digits); /* max number of decimal digits, fps   */

#define lrs_alloc_mp(a)
#define lrs_clear_mp(a)

lrs_mp_t lrs_alloc_mp_t(); /* dynamic allocation of lrs_mp                  */
lrs_mp_vector
lrs_alloc_mp_vector(long n); /* allocate lrs_mp_vector for n+1 lrs_mp numbers */
lrs_mp_matrix
lrs_alloc_mp_matrix(long m,
                    long n); /* allocate lrs_mp_matrix for m+1 x n+1 lrs_mp   */

void lrs_clear_mp_vector(lrs_mp_vector a, long n);
void lrs_clear_mp_matrix(lrs_mp_matrix a, long m, long n);

/*********************************************************/
/* Core library functions - depend on mp implementation  */
/******************************************************* */

long compare(lrs_mp a, lrs_mp b);     /* a ? b and returns -1,0,1 for <,=,> */
void gcd(lrs_mp u, lrs_mp v); /* returns u=gcd(u,v) destroying v         */
void mptodouble(lrs_mp a, double *x); /* convert lrs_mp to double */
long mptoi(lrs_mp a);                 /* convert lrs_mp to long integer */
char *mpgetstr10(char *, lrs_mp);     /* convert lrs_mp to string */
void reduce(lrs_mp Na, lrs_mp Da);    /* reduces Na Da by gcd(Na,Da) */

/*********************************************************/
/* Standard arithmetic & misc. functions                 */
/* should be independent of mp implementation            */
/******************************************************* */

long comprod(lrs_mp Na, lrs_mp Nb, lrs_mp Nc,
             lrs_mp Nd); /* +1 if Na*Nb > Nc*Nd,-1 if Na*Nb > Nc*Nd else 0 */
void getfactorial(lrs_mp factorial, long k); /* compute k factorial in lrs_mp */
void linrat(lrs_mp Na, lrs_mp Da, long ka, lrs_mp Nb, lrs_mp Db, long kb,
            lrs_mp Nc, lrs_mp Dc);
void lcm(lrs_mp a,
         lrs_mp b); /* a = least common multiple of a, b; b is saved  */
void notimpl(const char *s); /* bail out - help! */
void rattodouble(lrs_mp a, lrs_mp b,
                 double *x); /* convert lrs_mp rational to double          */
void reduceint(lrs_mp Na, lrs_mp Da); /* divide Na by Da and return it */
void reducearray(lrs_mp_vector p,
                 long n); /* find gcd of p[0]..p[n-1] and divide through by */
void scalerat(lrs_mp Na, lrs_mp Da, long ka); /* scales rational by ka */

/********************************/
/* Matrix and vector allocation */
/********************************/
void lrs_clear_mp_vector(lrs_mp_vector p, long n);
lrs_mp_vector lrs_alloc_mp_vector(long n);
void lrs_clear_mp_vector(lrs_mp_vector p, long n);
lrs_mp_matrix lrs_alloc_mp_matrix(long m, long n);
void lrs_clear_mp_matrix(lrs_mp_matrix p, long m, long n);
long lrs_mp_init(long dec_digits);

/**********************************/
/* Miscellaneous functions        */
/******************************** */

void lrs_getdigits(long *a, long *b); /* send digit information to user */

void stringcpy(char *s, char *t); /* copy t to s pointer version */

void *calloc();
void *malloc();
void *xcalloc(long n, long s, long l, const char *f);

void lrs_exit(int i);
void lrs_overflow(int i);
void lrsv2_overflow(int i);
/* end of  lrs_long.h (vertex enumeration using lexicographic reverse search) */
