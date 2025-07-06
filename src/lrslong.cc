/* lrslong.c      library code for lrs extended precision arithmetic */
/* Version 4.0, April 13, 2000                                       */
/* Copyright: David Avis 1999, avis@cs.mcgill.ca                     */

/* Derived from prs_single.c ( rational arithmetic for lrs and prs)  */
/* authored by  Ambros Marzetta    Revision 1.2  1998/05/27          */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// must be last or include errors
#include "lrslong.h"

__thread long lrs_digits;        /* max permitted no. of digits   */
__thread long lrs_record_digits; /* this is the biggest acheived so far.     */

#ifdef B128
__int128 MAXDm, MAXDl, MAXDa;
#endif

#define MAXINPUT 1000 /*max length of any input rational */

/* lrs_overflow routine should user supplied if not using lrslib.c */

/*
void lrs_overflow(int parm)
{
}
*/

void gcd(lrs_mp u, lrs_mp v)
/* Returns u=gcd(u,v) using classic Euclid's algorithm.
   v is destroyed.  Knuth, II, p.320 */

{
#ifndef B128
  unsigned long ul, vl, r;

  ul = labs(*u);
  vl = labs(*v);
#else
  __int128 ul, vl, r;
  ul = abs128(*u);
  vl = abs128(*v);
#endif
  if (ul == 0) {
    *u = vl;
    return;
  }

  while (vl != 0) {
    r = ul % vl;
    ul = vl;
    vl = r;
  }
  *u = ul;
} /* gcd */

void lcm(lrs_mp a,
         lrs_mp b) /* a = least common multiple of a, b; b is preserved */
{
  lrs_mp u, v;
  copy(u, a);
  copy(v, b);
  gcd(u, v);
  exactdivint(a, u, v); /* v=a/u   a contains remainder = 0 */
  mulint(v, b, a);
} /* end of lcm */

/***************************************************************/
/*                                                             */
/*     Package of routines for rational arithmetic             */
/*     (Built on top of package for multiprecision arithmetic  */
/*                                                             */
/***************************************************************/

void reduce(lrs_mp Na, lrs_mp Da) /* reduces Na/Da by gcd(Na,Da) */
{
  lrs_mp Nb, Db, Nc, Dc;
  copy(Nb, Na);
  copy(Db, Da);
  storesign(Nb, POS);
  storesign(Db, POS);
  copy(Nc, Na);
  copy(Dc, Da);
  gcd(Nb, Db); /* Nb is the gcd(Na,Da) */
  exactdivint(Nc, Nb, Na);
  exactdivint(Dc, Nb, Da);
}

void reduceint(lrs_mp Na, lrs_mp Da) /* divide Na by Da and return */
{
  lrs_mp Temp;
  copy(Temp, Na);
  exactdivint(Temp, Da, Na);
}

long comprod(lrs_mp Na, lrs_mp Nb, lrs_mp Nc,
             lrs_mp Nd) /* +1 if Na*Nb > Nc*Nd  */
                        /* -1 if Na*Nb < Nc*Nd  */
                        /*  0 if Na*Nb = Nc*Nd  */
{
  lrs_mp mc, md;
  itomp(ZERO, mc);
  itomp(ZERO, md); /*just to please some compilers */
  mulint(Na, Nb, mc);
  mulint(Nc, Nd, md);
  if (mp_greater(mc, md))
    return (1);
  if (mp_greater(md, mc))
    return (-1);
  return (0);
}

void linrat(lrs_mp Na, lrs_mp Da, long ka, lrs_mp Nb, lrs_mp Db, long kb,
            lrs_mp Nc, lrs_mp Dc)
/* computes Nc/Dc = ka*Na/Da  +kb* Nb/Db and reduces answer by gcd(Nc,Dc) */
{
  lrs_mp c;
  itomp(ZERO, c); /*just to please some compilers */
  mulint(Na, Db, Nc);
  mulint(Da, Nb, c);
  linint(Nc, ka, c, kb); /* Nc = (ka*Na*Db)+(kb*Da*Nb)  */
  mulint(Da, Db, Dc);    /* Dc =  Da*Db           */
  reduce(Nc, Dc);
}

/***************************************************************/
/*                                                             */
/*     Conversion and I/O functions                            */
/*                                                             */
/***************************************************************/

void mptodouble(lrs_mp a, double *x) /* convert lrs_mp to double */
{
  (*x) = (*a);
}

long mptoi(lrs_mp a) /* convert lrs_mp to long */
{
  return (*a);
}

/* return char * representation of a in base 10.
 * use out if non-NULL, otherwise allocate and return.
 */
char *mpgetstr10(char *out, lrs_mp a) {
  char *buf = NULL;
#ifndef B128
  int len = 0;
  len = snprintf(buf, 0, "%lld", *a);
  if (out != NULL)
    buf = out;
  else
    buf = static_cast<char *>(malloc(sizeof(char) * (len + 1)));
  sprintf(buf, "%lld", *a);
  return buf;
#else
  __int128 tmp = abs128(*a);
  int i, j = 0;
  int t[41];
  if (out != NULL)
    buf = out;
  else
    buf = static_cast<char *>(calloc(43, sizeof(char)));
  if (*a >= LLONG_MIN && *a <= LLONG_MAX) {
    sprintf(buf, "%lld", (long long)*a);
    return buf;
  }
  for (i = 0; tmp > 0; i++) {
    t[i] = tmp % 10;
    tmp = tmp / 10;
  }
  i--;
  if (*a < 0)
    buf[j++] = '-';
  while (i >= 0)
    buf[j++] = '0' + t[i--];
  return buf;
#endif
}

void rattodouble(lrs_mp a, lrs_mp b, double *x) /* convert lrs_mp rati
                                                   onal to double */

{
  double y;
  mptodouble(a, &y);
  mptodouble(b, x);
  *x = y / (*x);
}

/***************************************************************/
/*                                                             */
/*     Memory allocation functions                             */
/*                                                             */
/***************************************************************/
lrs_mp_t lrs_alloc_mp_t()
/* dynamic allocation of lrs_mp number */
{
  lrs_mp_t p;
  p = (lrs_mp_t)calloc(1, sizeof(lrs_mp));
  return p;
}

lrs_mp_vector lrs_alloc_mp_vector(long n)
/* allocate lrs_mp_vector for n+1 lrs_mp numbers */
{
  lrs_mp_vector p;
  long i;

  p = (lrs_mp_t *)CALLOC((n + 1), sizeof(lrs_mp_t));
  for (i = 0; i <= n; i++)
    p[i] = (lrs_mp_t)CALLOC(1, sizeof(lrs_mp));

  return p;
}

void lrs_clear_mp_vector(lrs_mp_vector p, long n)
/* free space allocated to p */
{
  long i;
  for (i = 0; i <= n; i++)
    free(p[i]);
  free(p);
  p = NULL;
}

lrs_mp_matrix lrs_alloc_mp_matrix(long m, long n)
/* allocate lrs_mp_matrix for m+1 x n+1 lrs_mp numbers */
{
  lrs_mp_matrix a;
  lrs_mp_t araw;
  int mp_width, row_width;
  int i, j;

  mp_width = lrs_digits + 1;
  row_width = (n + 1) * mp_width;

  araw = (lrs_mp_t)calloc((m + 1) * row_width, sizeof(lrs_mp));
  a = (lrs_mp_t **)calloc((m + 1), sizeof(lrs_mp_vector));

  for (i = 0; i < m + 1; i++) {
    a[i] = (lrs_mp_t *)calloc((n + 1), sizeof(lrs_mp_t));

    for (j = 0; j < n + 1; j++)
      a[i][j] = (araw + i * row_width + j * mp_width);
  }
  return a;
}

void lrs_clear_mp_matrix(lrs_mp_matrix p, long m, long n)
/* free space allocated to lrs_mp_matrix p */
{
  long i;

  /* p[0][0] is araw, the actual matrix storage address */

  free(p[0][0]);

  for (i = 0; i < m + 1; i++)
    free(p[i]);
  /* 2015.9.9 memory leak fix */
  free(p);
  p = NULL;
}

void lrs_getdigits(long *a, long *b) {
  /* send digit information to user */
  *a = DIG2DEC(ZERO);
  *b = DIG2DEC(ZERO);
  return;
}

void *xcalloc(long n, long s, long l, const char *f) {
  void *tmp;

  tmp = calloc(n, s);
  if (tmp == 0) {
    char buf[200];

    sprintf(buf, "\n\nFatal error on line %ld of %s", l, f);
    perror(buf);
    exit(1);
  }
  return tmp;
}

long lrs_mp_init(long dec_digits)
/* max number of decimal digits for the computation */
/* long int version                                 */
{

#ifdef B128
  MAXDl = 9223372036854775807L;       /* 2^63 - 1 */
  MAXDa = MAXDl * MAXDl;              /* should be 2^126 - 1 */
  MAXDm = 3611622602L;                /* sqrt(sqrt( 2^127 -1 ))         */
  MAXDm = MAXDm * MAXDm + 6000000000; /* too big to initialize directly */
#endif

  lrs_record_digits = 0;
  lrs_digits = 0; /* max permitted no. of digits   */
  return TRUE;
}

void notimpl(const char *s) {
  fflush(stdout);
  fprintf(stderr, "\nAbnormal Termination  %s\n", s);
  exit(1);
}

/***************************************************************/
/*                                                             */
/*     Misc. functions                                         */
/*                                                             */
/***************************************************************/

void reducearray(lrs_mp_vector p, long n)
/* find largest gcd of p[0]..p[n-1] and divide through */
{
  lrs_mp divisor;
  lrs_mp Temp;
  long i = 0L;

  while ((i < n) && zero(p[i]))
    i++;
  if (i == n)
    return;

  copy(divisor, p[i]);
  storesign(divisor, POS);
  i++;

  while (i < n) {
    if (!zero(p[i])) {
      copy(Temp, p[i]);
      storesign(Temp, POS);
      gcd(divisor, Temp);
    }
    i++;
  }

  /* reduce by divisor */
  for (i = 0; i < n; i++)
    if (!zero(p[i]))
      reduceint(p[i], divisor);
} /* end of reducearray */

void getfactorial(lrs_mp factorial, long k) /* compute k factorial
                                               in lrs_mp */
{
  lrs_mp temp;
  long i;
  itomp(ONE, factorial);
  for (i = 2; i <= k; i++) {
    itomp(i, temp);
    mulint(temp, factorial, factorial);
  }
} /* end of getfactorial */

void stringcpy(char *s, char *t) /*copy t to s pointer version */
{
  while (((*s++) = (*t++)) != '\0')
    ;
}
