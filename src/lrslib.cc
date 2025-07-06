/* lrslib.c     library code for lrs                     */

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

/*2022.1.22  lrs now has 3 modes controlled by these routines:
 fel_run      Fourier-Motzkin elimination, project or elimiminate options
*/

/* modified by Gary Roumanis and Skip Jordan for multithread compatability */
/* need to add a test for non-degenerate pivot step in reverse I guess?? */
/* Copyright: David Avis 2005,2022 avis@cs.mcgill.ca         */

#include <assert.h>
#include <exception>
#include <libgen.h>
#include <limits.h>
#include <setjmp.h>
#include <stdio.h>
#include <string.h>

#include "lrslib.h"

static __thread unsigned long dict_count, dict_limit, cache_tries, cache_misses;

/* Variables and functions global to this file only */

static __thread long lrs_checkpoint_seconds = 0;

static __thread long lrs_global_count = 0; /* Track how many lrs_dat records are
                                     allocated */

static __thread lrs_dat *lrs_global_list[MAX_LRS_GLOBALS + 1];

static lrs_dic *new_lrs_dic(long m, long d, long m_A);

static void cache_dict(lrs_dic **D_p, lrs_dat *global, long i, long j);
static long check_cache(lrs_dic **D_p, lrs_dat *global, long *i_p, long *j_p);

static void pushQ(lrs_dat *global, long m, long d, long m_A);

char *basename(char *path);

/*******************************/
/* functions  for external use */
/*******************************/

long lrs_getsolution(lrs_dic *P, lrs_dat *Q, lrs_mp_vector output, long col)
/* check if column indexed by col in this dictionary */
/* contains output                                   */
{

  long j; /* cobasic index     */

  lrs_mp_matrix A = P->A;
  long *Row = P->Row;

  if (col == ZERO) /* check for lexmin vertex */
    return lrs_getvertex(P, Q, output);

  /*  check for rays: negative in row 0 , positive if lponly */

  else if (!negative(A[0][col]))
    return FALSE;

  /*  and non-negative for all basic non decision variables */

  j = Q->lastdv + 1;
  while (j <= P->m && !negative(A[Row[j]][col]))
    j++;

  if (j <= P->m)
    return FALSE;
  /* 2021.5.21 check for max depth output */
  if (lexmin(P, Q, col)) {
    if (P->depth > Q->count[8])
      Q->count[8] = P->depth;
    return lrs_getray(P, Q, col, Q->n, output);
  }

  return FALSE; /* no more output in this dictionary */

} /* end of lrs_getsolution */

long lrs_init() /* returns TRUE if successful, else FALSE */
{
  if (!lrs_mp_init(0)) /* initialize arithmetic */
    return FALSE;
  lrs_global_count = 0;
  lrs_checkpoint_seconds = 0;
  return TRUE;
}

/***********************************/
/* allocate and initialize lrs_dat */
/***********************************/
lrs_dat *lrs_alloc_dat() {
  lrs_dat *Q;
  long i;

  if (lrs_global_count >= MAX_LRS_GLOBALS) {
    return NULL;
  }

  Q = (lrs_dat *)malloc(sizeof(lrs_dat));
  if (Q == NULL)
    return Q; /* failure to allocate */

  lrs_global_list[lrs_global_count] = Q;
  lrs_global_count++;

  Q->m = 0L;
  Q->n = 0L;
  Q->nlinearity = 0L;
  Q->nredundcol = 0L;
  for (i = 0; i < 10; i++) {
    Q->count[i] = 0L;
  }
  Q->count[2] = 1L; /* basis counter */
                    /* initialize flags */

  Q->redundcol = NULL;
  Q->inequality = NULL;
  Q->linearity = NULL;
  Q->minratio = NULL;
  Q->temparray = NULL;
  Q->olddic = NULL;

  lrs_alloc_mp(Q->Nvolume);
  lrs_alloc_mp(Q->Dvolume);
  itomp(ZERO, Q->Nvolume);
  itomp(ONE, Q->Dvolume);

  return Q;
} /* end of allocate and initialize lrs_dat */

/* In lrs_getfirstbasis and lrs_getnextbasis we use D instead of P */
/* since the dictionary P may change, ie. &P in calling routine    */

#define D (*D_p)

long lrs_getfirstbasis(lrs_dic **D_p, lrs_dat *Q, lrs_mp_matrix *Lin,
                       long no_output)
/* gets first basis, FALSE if none              */
/* P may get changed if lin. space Lin found    */
/* no_output is TRUE supresses output headers   */
/* 2017.12.22  could use no_output=2 to get early exit for criss-cross method */
{
  long i, j, k;

  /* assign local variables to structures */

  lrs_mp_matrix A;
  long *B, *C, *Col;
  long *inequality;
  long *linearity;
  long m, d, lastdv, nlinearity, nredundcol;

  m = D->m;
  d = D->d;

  lastdv = Q->lastdv;

  nredundcol = 0L; /* will be set after getabasis        */
  nlinearity =
      Q->nlinearity; /* may be reset if new linearity read or in getabasis*/
  linearity = Q->linearity;

  A = D->A;
  B = D->B;
  C = D->C;
  Col = D->Col;
  inequality = Q->inequality;

  /**2020.6.17 with extract just select columns and quit */

  /* default is to look for starting cobasis using linearies first, then     */
  /* filling in from last rows of input as necessary                         */
  /* linearity array is assumed sorted here                                  */
  /* note if restart/given start inequality indices already in place         */
  /* from nlinearity..d-1                                                    */
  for (i = 0; i < nlinearity; i++) /* put linearities first in the order */
    inequality[i] = linearity[i];

  k = nlinearity;

  for (i = m; i >= 1; i--) {
    j = 0;
    while (j < k && inequality[j] != i)
      j++; /* see if i is in inequality  */
    if (j == k)
      inequality[k++] = i;
  }

  for (j = 0; j <= d; j++)
    itomp(ZERO, A[0][j]);

  /* Now we pivot to standard form, and then find a primal feasible basis */
  /* Note these steps MUST be done, even if restarting, in order to get */
  /* the same index/inequality correspondance we had for the original prob. */
  /* The inequality array is used to give the insertion order */

  if (!getabasis(D, Q, inequality))
    return FALSE;
  /* bug fix 2009.12.2 */
  nlinearity =
      Q->nlinearity; /*may have been reset if some lins are redundant*/

  /* 2020.2.2 */
  /* extract option asked to remove all linearities and output the reduced A
   * matrix */
  /* should be followed by redund to get minimum representation */

  nredundcol = Q->nredundcol;
  lastdv = Q->lastdv;
  d = D->d;

  /* Reset up the inequality array to remember which index is which input
   * inequality */
  /* inequality[B[i]-lastdv] is row number of the inequality with index B[i] */
  /* inequality[C[i]-lastdv] is row number of the inequality with index C[i] */

  for (i = 1; i <= m; i++)
    inequality[i] = i;
  if (nlinearity > 0) /* some cobasic indices will be removed */
  {
    for (i = 0; i < nlinearity; i++) /* remove input linearity indices */
      inequality[linearity[i]] = 0;
    k = 1; /* counter for linearities         */
    for (i = 1; i <= m - nlinearity; i++) {
      while (k <= m && inequality[k] == 0)
        k++; /* skip zeroes in corr. to linearity */
      inequality[i] = inequality[k++];
    }
  }
  /* end if linearity */

  if (nredundcol > 0) {
    const unsigned int Qn = Q->n;
    *Lin = lrs_alloc_mp_matrix(nredundcol, Qn);

    for (i = 0; i < nredundcol; i++) {

      lrs_getray(D, Q, Col[0], D->C[0] + i,
                  (*Lin)[i]); /* adjust index for deletions */

      if (!removecobasicindex(D, Q, 0L)) {
        lrs_clear_mp_matrix(*Lin, nredundcol, Qn);
        return FALSE;
      }
    }

  } /* end if nredundcol > 0 */

  /*2017.12.22   If you want to do criss-cross now is the time ! */

  /* Do dual pivots to get primal feasibility */
  if (!primalfeasible(D, Q)) {
    return FALSE;
  }

  /* Now solve LP if objective function was given */

  for (j = 1; j <= d; j++) {
    copy(A[0][j], D->det);
    storesign(A[0][j], NEG);
  }

  itomp(ZERO, A[0][0]); /* zero optimum objective value */

  /* reindex basis to 0..m if necessary */
  /* we use the fact that cobases are sorted by index value */

  while (C[0] <= m) {
    i = C[0];
    j = inequality[B[i] - lastdv];
    inequality[B[i] - lastdv] = inequality[C[0] - lastdv];
    inequality[C[0] - lastdv] = j;
    C[0] = B[i];
    B[i] = i;
    reorder1(C, Col, ZERO, d);
  }

  /* Check to see if necessary to resize */
  /* bug fix 2018.6.7 new value of d required below */
  if (Q->n - 1 > D->d)
    *D_p = resize(D, Q);
  return TRUE;
}
/********* end of lrs_getfirstbasis  ***************/

/*****************************************/
/* getnextbasis in reverse search order  */
/*****************************************/

long lrs_getnextbasis(lrs_dic **D_p, lrs_dat *Q)
{
  /* assign local variables to structures */
  long i = 0L, j = 0L;
  long m = D->m;
  long d = D->d;

  long backtrack = FALSE;

  while ((j < d) || (D->B[m] != m)) /*main while loop for getnextbasis */
  {
    if (D->depth >= MAXD) {
      backtrack = TRUE;
    }

    if (backtrack) /* go back to prev. dictionary, restore i,j */
    {
      backtrack = FALSE;

      if (check_cache(D_p, Q, &i, &j)) {
      } else {
        D->depth--;
        selectpivot(D, Q, &i, &j);
        pivot(D, Q, i, j);
        update(D, Q, &i, &j); /*Update B,C,i,j */
      }

      j++; /* go to next column */
    } /* end of if backtrack  */

    if (D->depth < -MAXD)
      break;

    /* try to go down tree */

    /* 2011.7.14 patch */
    while ((j < d) && (!reverse(D, Q, &i, j)))
      j++;
    if (j == d)
      backtrack = TRUE;
    else
    /*reverse pivot found */
    {
      cache_dict(D_p, Q, i, j);
      /* Note that the next two lines must come _after_ the
         call to cache_dict */

      D->depth++;

      pivot(D, Q, i, j);
      update(D, Q, &i, &j); /*Update B,C,i,j */

      D->lexflag = lexmin(D, Q, ZERO); /* see if lexmin basis */
      Q->count[2]++;

      return TRUE; /*return new dictionary */
    }

  } /* end of main while loop for getnextbasis */
  return FALSE; /* done, no more bases */
} /*end of lrs_getnextbasis */

/*************************************/
/* print out one line of output file */
/*************************************/
long lrs_getvertex(lrs_dic *P, lrs_dat *Q, lrs_mp_vector output)
/*Print out current vertex if it is lexmin and return it in output */
/* return FALSE if no output generated  */
{
  long i;
  long ind;  /* output index                                  */
  long ired; /* counts number of redundant columns            */
             /* assign local variables to structures */
  long *redundcol = Q->redundcol;

  long lexflag;

  lexflag = P->lexflag;
  if (lexflag) {
    ++(Q->count[1]);
    /* 2021.5.21 check for max depth output */
    if (P->depth > Q->count[8])
      Q->count[8] = P->depth;
  }

  if (!lexflag) /* not lexmin, and not printing forced */
    return FALSE;

  /* copy column 0 to output */

  i = 1;
  ired = 0;
  copy(output[0], P->det);

  for (ind = 1; ind < Q->n; ind++) /* extract solution */

    if ((ired < Q->nredundcol) && (redundcol[ired] == ind))
    /* column was deleted as redundant */
    {
      itomp(ZERO, output[ind]);
      ired++;
    } else
    /* column not deleted as redundant */
    {
      getnextoutput(P, Q, i, ZERO, output[ind]);
      i++;
    }

  reducearray(output, Q->n);
  if (lexflag && one(output[0]))
    ++Q->count[4]; /* integer vertex */

  return TRUE;
} /* end of lrs_getvertex */

long lrs_getray(lrs_dic *P, lrs_dat *Q, long col, long redcol,
                lrs_mp_vector output)
/*Print out solution in col and return it in output   */
/* return FALSE if no output generated in column col  */
{
  long i;
  long ind;  /* output index                                  */
  long ired; /* counts number of redundant columns            */
             /* assign local variables to structures */
  long *redundcol = Q->redundcol;
  long *count = Q->count;
  long n = Q->n;

  if (redcol == n) {
    ++count[0];
  }

  i = 1;
  ired = 0;

  for (ind = 0; ind < n; ind++) /* print solution */
  {
    if (ind == 0) /* must have a ray, set first column to zero */
      itomp(ZERO, output[0]);

    else if ((ired < Q->nredundcol) && (redundcol[ired] == ind))
    /* column was deleted as redundant */
    {
      if (redcol == ind) /* true for linearity on this cobasic index */
        /* we print reduced determinant instead of zero */
        copy(output[ind], P->det);
      else
        itomp(ZERO, output[ind]);
      ired++;
    } else
    /* column not deleted as redundant */
    {
      getnextoutput(P, Q, i, col, output[ind]);
      i++;
    }
  }
  reducearray(output, n);

  return TRUE;
} /* end of lrs_getray */

void getnextoutput(lrs_dic *P, lrs_dat *Q, long i, long col, lrs_mp out)
/* get A[B[i][col] and copy to out */
{
  long row;
  lrs_mp_matrix A = P->A;
  long *Row = P->Row;

  row = Row[i];

  copy(out, A[row][col]);

} /* end of getnextoutput */

/*********************************/
/* Internal functions            */
/*********************************/
/* Basic Dictionary functions    */
/******************************* */

long reverse(lrs_dic *P, lrs_dat *Q, long *r, long s)
/*  find reverse indices  */
/* TRUE if B[*r] C[s] is a reverse lexicographic pivot */
{
  long i, j, row, col;

  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long d = P->d;

  col = Col[s];
  if (!negative(A[0][col])) {
    Q->minratio[P->m] = 0; /* 2011.7.14 */
    return (FALSE);
  }

  *r = lrs_ratio(P, Q, col);
  if (*r == 0) /* we have a ray */
  {
    Q->minratio[P->m] = 0; /* 2011.7.14 */
    return (FALSE);
  }

  row = Row[*r];

  /* check cost row after "pivot" for smaller leaving index    */
  /* ie. j s.t.  A[0][j]*A[row][col] < A[0][col]*A[row][j]     */
  /* note both A[row][col] and A[0][col] are negative          */

  for (i = 0; i < d && C[i] < B[*r]; i++)
    if (i != s) {
      j = Col[i];
      if (positive(A[0][j]) ||
          negative(A[row][j])) /*or else sign test fails trivially */
        if ((!negative(A[0][j]) && !positive(A[row][j])) ||
            comprod(A[0][j], A[row][col], A[0][col], A[row][j]) ==
                -1) {            /*+ve cost found */
          Q->minratio[P->m] = 0; /* 2011.7.14 */

          return (FALSE);
        }
    }
  return (TRUE);
} /* end of reverse */

long selectpivot(lrs_dic *P, lrs_dat *Q, long *r, long *s)
/* select pivot indices using lexicographic rule   */
/* returns TRUE if pivot found else FALSE          */
/* pivot variables are B[*r] C[*s] in locations Row[*r] Col[*s] */
{
  long j, col;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *Col = P->Col;
  long d = P->d;

  *r = 0;
  *s = d;
  j = 0;

  /*find positive cost coef */
  while ((j < d) && (!positive(A[0][Col[j]])))
    j++;

  if (j < d) /* pivot column found! */
  {
    *s = j;
    col = Col[j];

    /*find min index ratio */
    *r = lrs_ratio(P, Q, col);
    if (*r != 0)
      return (TRUE); /* unbounded if *r=0 */
  }
  return (FALSE);
} /* end of selectpivot        */
/******************************************************* */

void pivot(lrs_dic *P, lrs_dat *Q, long bas, long cob)
/* Qpivot routine for array A              */
/* indices bas, cob are for Basis B and CoBasis C    */
/* corresponding to row Row[bas] and column       */
/* Col[cob]   respectively                       */
{
  long r, s;
  long i, j;

  /* assign local variables to structures */

  lrs_mp_matrix A = P->A;
  long *Row = P->Row;
  long *Col = P->Col;
  long d, m_A;

#ifndef LRSLONG
  lrs_mp Ns, Nt;
  lrs_alloc_mp(Ns);
  lrs_alloc_mp(Nt);
#endif
  lrs_mp Ars;
  lrs_alloc_mp(Ars);

  d = P->d;
  m_A = P->m_A;
  Q->count[3]++; /* count the pivot */

  r = Row[bas];
  s = Col[cob];

  /* Ars=A[r][s]    */
  copy(Ars, A[r][s]);
  storesign(P->det, sign(Ars)); /*adjust determinant to new sign */

  for (i = 0; i <= m_A; i++)
    if (i != r)
      for (j = 0; j <= d; j++)
        if (j != s) {
          /*          A[i][j]=(A[i][j]*Ars-A[i][s]*A[r][j])/P->det; */

#ifdef LRSLONG
          qpiv(A[i][j], Ars, A[i][s], A[r][j], P->det);
#else
          mulint(A[i][j], Ars, Nt);
          mulint(A[i][s], A[r][j], Ns);
          decint(Nt, Ns);
          exactdivint(Nt, P->det, A[i][j]);
#endif
        } /* end if j ....  */

  if (sign(Ars) == POS) {
    for (j = 0; j <= d; j++) /* no need to change sign if Ars neg */
      /*   A[r][j]=-A[r][j];              */
      if (!zero(A[r][j]))
        changesign(A[r][j]);
  } /* watch out for above "if" when removing this "}" ! */
  else
    for (i = 0; i <= m_A; i++)
      if (!zero(A[i][s]))
        changesign(A[i][s]);

  /*  A[r][s]=P->det;                  */

  copy(A[r][s], P->det); /* restore old determinant */
  copy(P->det, Ars);
  storesign(P->det, POS); /* always keep positive determinant */

  /* set the new rescaled objective function value */

  mulint(P->det, Q->Lcm[0], P->objden);
  mulint(Q->Gcd[0], A[0][0], P->objnum);

  changesign(P->objnum);

  if (zero(P->objnum))
    storesign(P->objnum, POS);
  else
    reduce(P->objnum, P->objden);

#ifndef LRSLONG
  lrs_clear_mp(Ns);
  lrs_clear_mp(Nt);
#endif
  lrs_clear_mp(Ars);
} /* end of pivot */

long primalfeasible(lrs_dic *P, lrs_dat *Q)
/* Do dual pivots to get primal feasibility */
/* Note that cost row is all zero, so no ratio test needed for Dual Bland's rule
 */
{
  long primalinfeasible = TRUE;
  long i, j;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *Row = P->Row;
  long *Col = P->Col;
  long m, d, lastdv;
  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;

  /*temporary: try to get new start after linearity */

  while (primalinfeasible) {
    i = lastdv + 1;
    while (i <= m && !negative(A[Row[i]][0]))
      i++;
    if (i <= m) {
      j = 0; /*find a positive entry for in row */
      while (j < d && !positive(A[Row[i]][Col[j]]))
        j++;
      if (j >= d)
        return (FALSE); /* no positive entry */
      pivot(P, Q, i, j);
      update(P, Q, &i, &j);
    } else
      primalinfeasible = FALSE;
  } /* end of while primalinfeasibile */
  return (TRUE);
} /* end of primalfeasible */

long lrs_solvelp(lrs_dic *P, lrs_dat *Q, long maximize)
/* Solve primal feasible lp by Dantzig`s rule and lexicographic ratio test */
/* return TRUE if bounded, FALSE if unbounded                              */
{
  long i = 0L;
  /* assign local variables to structures */
  long d = P->d;

  if (0 < d && i == 0) /* selectpivot gives information on unbounded solution */
  {
    return FALSE;
  }
  return TRUE;
} /* end of lrs_solvelp  */

long getabasis(lrs_dic *P, lrs_dat *Q, long order[])

/* Pivot Ax<=b to standard form */
/*Try to find a starting basis by pivoting in the variables x[1]..x[d]        */
/*If there are any input linearities, these appear first in order[]           */
/* Steps: (a) Try to pivot out basic variables using order                    */
/*            Stop if some linearity cannot be made to leave basis            */
/*        (b) Permanently remove the cobasic indices of linearities           */
/*        (c) If some decision variable cobasic, it is a linearity,           */
/*            and will be removed.                                            */

{
  long i, j, k;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long *linearity = Q->linearity;
  long *redundcol = Q->redundcol;
  long m, d, nlinearity;
  long nredundcol = 0L; /* will be calculated here */
  m = P->m;
  d = P->d;
  nlinearity = Q->nlinearity;

  for (j = 0l; j < m; j++) {
    i = 0l;
    while (i <= m && B[i] != d + order[j])
      i++;                       /* find leaving basis index i */
    if (j < nlinearity && i > m) /* cannot pivot linearity to cobasis */
    {
      return FALSE;
    }
    if (i <= m) { /* try to do a pivot */
      k = 0l;
      while (C[k] <= d && zero(A[Row[i]][Col[k]])) {
        k++;
      }
      if (C[k] <= d) {

        pivot(P, Q, i, k);
        update(P, Q, &i, &k);
      } else if (j < nlinearity) { /* cannot pivot linearity to cobasis */
        if (zero(A[Row[i]][0])) {
          linearity[j] = 0l;
        } else {
          return FALSE;
        }
      }
    }
  }

  /* update linearity array to get rid of redundancies */
  i = 0;
  k = 0; /* counters for linearities         */
  while (k < nlinearity) {
    while (k < nlinearity && linearity[k] == 0)
      k++;
    if (k < nlinearity)
      linearity[i++] = linearity[k++];
  }

  nlinearity = i;
  /* bug fix, 2009.6.27 */ Q->nlinearity = i;

  /* column dependencies now can be recorded  */
  /* redundcol contains input column number 0..n-1 where redundancy is */
  k = 0;
  while (k < d && C[k] <= d) {
    if (C[k] <= d) { /* decision variable still in cobasis */
      redundcol[nredundcol++] = C[k];
    }
    k++;
  }

  /* now we know how many decision variables remain in problem */
  Q->nredundcol = nredundcol;
  Q->lastdv = d - nredundcol;

  /* Remove linearities from cobasis for rest of computation */
  /* This is done in order so indexing is not screwed up */

  for (i = 0; i < nlinearity; i++) { /* find cobasic index */
    k = 0;
    while (k < d && C[k] != linearity[i] + d)
      k++;
    if (k >= d) {
      return FALSE;
    }
    if (!removecobasicindex(P, Q, k))
      return FALSE;
    d = P->d;
  }

  return TRUE;
} /*  end of getabasis */

long removecobasicindex(lrs_dic *P, lrs_dat *Q, long k)
/* remove the variable C[k] from the problem */
/* used after detecting column dependency    */
{
  long i, j, cindex, deloc;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Col = P->Col;
  long m, d;
  m = P->m;
  d = P->d;

  cindex = C[k];  /* cobasic index to remove              */
  deloc = Col[k]; /* matrix column location to remove     */

  for (i = 1; i <= m; i++) /* reduce basic indices by 1 after index */
    if (B[i] > cindex)
      B[i]--;

  for (j = k; j < d; j++) /* move down other cobasic variables    */
  {
    C[j] = C[j + 1] - 1; /* cobasic index reduced by 1           */
    Col[j] = Col[j + 1];
  }

  if (deloc != d) {
    /* copy col d to deloc */
    for (i = 0; i <= m; i++)
      copy(A[i][deloc], A[i][d]);

    /* reassign location for moved column */
    j = 0;
    while (Col[j] != d)
      j++;

    Col[j] = deloc;
  }

  P->d--;
  return TRUE;
} /* end of removecobasicindex */

lrs_dic *resize(lrs_dic *P, lrs_dat *Q)
{
  lrs_dic *P1; /* to hold new dictionary in case of resizing */

  long i, j;
  long m, d, m_A;

  m = P->m;
  d = P->d;
  m_A = P->m_A;

  /* get new dictionary record */

  P1 = new_lrs_dic(m, d, m_A);

  /* copy data from P to P1    */
  P1->i = P->i;
  P1->j = P->j;
  P1->depth = P->depth;
  P1->m = P->m;
  P1->d = P1->d_orig = d;
  P1->lexflag = P->lexflag;
  P1->m_A = P->m_A;
  copy(P1->det, P->det);
  copy(P1->objnum, P->objnum);
  copy(P1->objden, P->objden);

  for (i = 0; i <= m; i++) {
    P1->B[i] = P->B[i];
    P1->Row[i] = P->Row[i];
  }
  for (i = 0; i <= m_A; i++) {
    for (j = 0; j <= d; j++)
      copy(P1->A[i][j], P->A[i][j]);
  }

  for (j = 0; j <= d; j++) {
    P1->Col[j] = P->Col[j];
    P1->C[j] = P->C[j];
  }

  lrs_free_dic(P, Q);

  /* Reassign cache pointers */

  Q->Qhead = P1;
  Q->Qtail = P1;
  P1->next = P1;
  P1->prev = P1;

  return P1;
}
/********* resize                    ***************/

long lrs_ratio(lrs_dic *P, lrs_dat *Q, long col) /*find lex min. ratio */
/* find min index ratio -aig/ais, ais<0 */
/* if multiple, checks successive basis columns */
/* recoded Dec 1997                     */
{
  long i, j, comp, ratiocol, basicindex, start, nstart, cindex, bindex;
  long firstime; /*For ratio test, true on first pass,else false */
  lrs_mp Nmin, Dmin;
  long degencount, ndegencount;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *Row = P->Row;
  long *Col = P->Col;
  long *minratio = Q->minratio;
  long m, d, lastdv;

  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;

  nstart = 0;
  ndegencount = 0;
  degencount = 0;
  minratio[P->m] = 1; /*2011.7.14 non-degenerate pivot flag */

  for (j = lastdv + 1; j <= m; j++) {
    /* search rows with negative coefficient in dictionary */
    /*  minratio contains indices of min ratio cols        */
    if (negative(A[Row[j]][col])) {
      minratio[degencount++] = j;
      if (zero(A[Row[j]][0]))
        minratio[P->m] = 0; /*2011.7.14 degenerate pivot flag */
    }
  } /* end of for loop */
  if (degencount == 0)
    return (degencount); /* non-negative pivot column */

  lrs_alloc_mp(Nmin);
  lrs_alloc_mp(Dmin);
  ratiocol = 0;   /* column being checked, initially rhs */
  start = 0;      /* starting location in minratio array */
  bindex = d + 1; /* index of next basic variable to consider */
  cindex = 0;     /* index of next cobasic variable to consider */
  basicindex =
      d; /* index of basis inverse for current ratio test, except d=rhs test */
  while (degencount > 1) /*keep going until unique min ratio found */
  {
    if (B[bindex] == basicindex) /* identity col in basis inverse */
    {
      if (minratio[start] == bindex)
      /* remove this index, all others stay */
      {
        start++;
        degencount--;
      }
      bindex++;
    } else
    /* perform ratio test on rhs or column of basis inverse */
    {
      firstime = TRUE;
      /*get next ratio column and increment cindex */
      if (basicindex != d)
        ratiocol = Col[cindex++];
      for (j = start; j < start + degencount; j++) {
        i = Row[minratio[j]]; /* i is the row location of the next basic
                                 variable */
        comp = 1;             /* 1:  lhs>rhs;  0:lhs=rhs; -1: lhs<rhs */
        if (firstime)
          firstime = FALSE; /*force new min ratio on first time */
        else {
          if (positive(Nmin) || negative(A[i][ratiocol])) {
            if (negative(Nmin) || positive(A[i][ratiocol]))
              comp = comprod(Nmin, A[i][col], A[i][ratiocol], Dmin);
            else
              comp = -1;
          }

          else if (zero(Nmin) && zero(A[i][ratiocol]))
            comp = 0;

          if (ratiocol == ZERO)
            comp = -comp; /* all signs reversed for rhs */
        }
        if (comp == 1) { /*new minimum ratio */
          nstart = j;
          copy(Nmin, A[i][ratiocol]);
          copy(Dmin, A[i][col]);
          ndegencount = 1;
        } else if (comp == 0) /* repeated minimum */
          minratio[nstart + ndegencount++] = minratio[j];

      } /* end of  for (j=start.... */
      degencount = ndegencount;
      start = nstart;
    } /* end of else perform ratio test statement */
    basicindex++; /* increment column of basis inverse to check next */
  } /*end of while loop */
  lrs_clear_mp(Nmin);
  lrs_clear_mp(Dmin);
  return (minratio[start]);
} /* end of ratio */

long lexmin(lrs_dic *P, lrs_dat *Q, long col)
/*test if basis is lex-min for vertex or ray, if so TRUE */
/* FALSE if a_r,g=0, a_rs !=0, r > s          */
{
  /*do lexmin test for vertex if col=0, otherwise for ray */
  long r, s, i, j;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long m = P->m;
  long d = P->d;

  for (i = Q->lastdv + 1; i <= m; i++) {
    r = Row[i];
    if (zero(A[r][col])) /* necessary for lexmin to fail */
      for (j = 0; j < d; j++) {
        s = Col[j];
        if (B[i] > C[j]) /* possible pivot to reduce basis */
        {
          if (zero(A[r][0])) /* no need for ratio test, any pivot feasible */
          {
            if (!zero(A[r][s]))
              return (FALSE);
          } else if (negative(A[r][s]) && ismin(P, Q, r, s)) {
            return (FALSE);
          }
        } /* end of if B[i] ... */
      }
  }
  return (TRUE);
} /* end of lexmin */

long ismin(lrs_dic *P, lrs_dat *Q, long r, long s)
/*test if A[r][s] is a min ratio for col s */
{
  long i;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long m_A = P->m_A;

  for (i = 1; i <= m_A; i++)
    if ((i != r) && negative(A[i][s]) &&
        comprod(A[i][0], A[r][s], A[i][s], A[r][0])) {
      return (FALSE);
    }

  return (TRUE);
}

void update(lrs_dic *P, lrs_dat *Q, long *i, long *j)
/*update the B,C arrays after a pivot */
/*   involving B[bas] and C[cob]           */
{

  long leave, enter;
  /* assign local variables to structures */
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long m = P->m;
  long d = P->d;

  leave = B[*i];
  enter = C[*j];
  B[*i] = enter;
  reorder1(B, Row, *i, m + 1);
  C[*j] = leave;
  reorder1(C, Col, *j, d);
  /* restore i and j to new positions in basis */
  for (*i = 1; B[*i] != enter; (*i)++)
    ; /*Find basis index */
  for (*j = 0; C[*j] != leave; (*j)++)
    ; /*Find co-basis index */
} /* end of update */

/*********************************************************/
/*                 Miscellaneous                         */
/******************************************************* */

void reorder(long a[], long range)
/*reorder array in increasing order with one misplaced element */
{
  long i, temp;
  for (i = 0; i < range - 1; i++)
    if (a[i] > a[i + 1]) {
      temp = a[i];
      a[i] = a[i + 1];
      a[i + 1] = temp;
    }
  for (i = range - 2; i >= 0; i--)
    if (a[i] > a[i + 1]) {
      temp = a[i];
      a[i] = a[i + 1];
      a[i + 1] = temp;
    }

} /* end of reorder */

void reorder1(long a[], long b[], long newone, long range)
/*reorder array a in increasing order with one misplaced element at index newone
 */
/*elements of array b are updated to stay aligned with a */
{
  long temp;
  while (newone > 0 && a[newone] < a[newone - 1]) {
    temp = a[newone];
    a[newone] = a[newone - 1];
    a[newone - 1] = temp;
    temp = b[newone];
    b[newone] = b[newone - 1];
    b[--newone] = temp;
  }
  while (newone < range - 1 && a[newone] > a[newone + 1]) {
    temp = a[newone];
    a[newone] = a[newone + 1];
    a[newone + 1] = temp;
    temp = b[newone];
    b[newone] = b[newone + 1];
    b[++newone] = temp;
  }
} /* end of reorder1 */

/***************************************************/
/* Routines for redundancy checking                */
/***************************************************/

long checkredund(lrs_dic *P, lrs_dat *Q)
/* Solve primal feasible lp by least subscript and lex min basis method */
/* to check redundancy of a row in objective function                   */
/* return 0=nonredundant -1=strict redundant(interior) 1=non-strict redundant*/
{
  lrs_mp Ns, Nt;
  long i, j;
  long r, s;

  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *Row, *Col;
  long d = P->d;

  lrs_alloc_mp(Ns);
  lrs_alloc_mp(Nt);
  Row = P->Row;
  Col = P->Col;
  while (selectpivot(P, Q, &i, &j)) {
    Q->count[2]++;

    /* sign of new value of A[0][0]            */
    /* is      A[0][s]*A[r][0]-A[0][0]*A[r][s] */

    r = Row[i];
    s = Col[j];

    mulint(A[0][s], A[r][0], Ns);
    mulint(A[0][0], A[r][s], Nt);

    if (mp_greater(Ns, Nt)) {
      lrs_clear_mp(Ns);
      lrs_clear_mp(Nt);
      return 0; /* non-redundant */
    }

    pivot(P, Q, i, j);
    update(P, Q, &i, &j); /*Update B,C,i,j */
  }
  lrs_clear_mp(Ns);
  lrs_clear_mp(Nt);

  if (j < d && i == 0) /* unbounded is also non-redundant */
    return 0;

  /* 2020.6.8 check for strict redundancy and return -1 if so */

  if (negative(P->A[0][0]))
    return -1;
  else
    return 1;

} /* end of checkredund  */

long checkcobasic(lrs_dic *P, lrs_dat *Q, long index)
/* TRUE if index is cobasic and nonredundant                         */
/* FALSE if basic, or degen. cobasic, where it will get pivoted out  */

{

  /* assign local variables to structures */

  lrs_mp_matrix A = P->A;
  long *C, *Row, *Col;
  long d = P->d;
  long m = P->m;
  long i = 0;
  long j = 0;
  long s;

  C = P->C;
  Row = P->Row;
  Col = P->Col;

  while ((j < d) && C[j] != index)
    j++;

  if (j == d)
    return FALSE; /* not cobasic index */

  /* index is cobasic */

  s = Col[j];
  i = Q->lastdv + 1;

  while ((i <= m) && (zero(A[Row[i]][s]) || !zero(A[Row[i]][0])))
    i++;

  if (i > m) {
    return TRUE;
  }

  pivot(P, Q, i, j);
  update(P, Q, &i, &j); /*Update B,C,i,j */

  return FALSE; /*index is no longer cobasic */

} /* end of checkcobasic */

long checkindex(lrs_dic *P, lrs_dat *Q, long index)
/* 0 if index is non-redundant inequality    */
/*-1 if index is strict redundant inequality */
/* 1 if index is non-strict redundant ine    */
/* 2 if index is input linearity             */
/*NOTE: row is returned all zero if redundant!! */
{
  long i, j;

  lrs_mp_matrix A = P->A;
  long *Row = P->Row;
  long *B = P->B;
  long d = P->d;
  long m = P->m;
  long zeroonly = 0;

  if (index < 0)
  {
    zeroonly = 1;
    index = -index;
  }

  /* each slack index must be checked for redundancy */
  /* if in cobasis, it is pivoted out if degenerate */
  /* else it is non-redundant                       */

  if (checkcobasic(P, Q, index)) {
    return ZERO;
  }
  /* index is basic   */
  j = 1;
  while ((j <= m) && (B[j] != index))
    j++;

  i = Row[j];

  /* copy row i to cost row, and set it to zero */

  for (j = 0; j <= d; j++) {
    copy(A[0][j], A[i][j]);
    changesign(A[0][j]);
    itomp(ZERO, A[i][j]);
  }
  if (zeroonly)
    return 1;

  /*2020.6.6 new test for strict redundancy */

  j = checkredund(P, Q);
  if (j != 0)
    return j;

  /* non-redundant, copy back and change sign */

  for (j = 0; j <= d; j++) {
    copy(A[i][j], A[0][j]);
    changesign(A[i][j]);
  }

  return 0;

} /* end of checkindex */

/***************************************************************/
/*                                                             */
/*     Routines for caching, allocating etc.                   */
/*                                                             */
/***************************************************************/

/* From here mostly Bremner's handiwork */

static void cache_dict(lrs_dic **D_p, lrs_dat *global, long i, long j) {

  if (dict_limit > 1) {
    /* save row, column indicies */
    (*D_p)->i = i;
    (*D_p)->j = j;

    /* Make a new, blank spot at the end of the queue to copy into */

    pushQ(global, (*D_p)->m, (*D_p)->d, (*D_p)->m_A);

    /*2019.6.7 This ought not to happen but it does */

    if (global->Qtail == *D_p)
      return;

    copy_dict(global, global->Qtail, *D_p); /* Copy current dictionary */
  }
  *D_p = global->Qtail;
}

void copy_dict(lrs_dat *global, lrs_dic *dest, lrs_dic *src) {
  long m = src->m;
  long m_A = src->m_A; /* number of rows in A */
  long d = src->d;
  long r, s;

  assert(dest != src);

  for (r = 0; r <= m_A; r++)
    for (s = 0; s <= d; s++)
      copy(dest->A[r][s], src->A[r][s]);

  dest->i = src->i;
  dest->j = src->j;
  dest->m = m;
  dest->d = d;
  dest->d_orig = src->d_orig;
  dest->m_A = src->m_A;

  dest->depth = src->depth;
  dest->lexflag = src->lexflag;

  copy(dest->det, src->det);
  copy(dest->objnum, src->objnum);
  copy(dest->objden, src->objden);

  memcpy(dest->B, src->B, (m + 1) * sizeof(long));
  memcpy(dest->C, src->C, (d + 1) * sizeof(long));
  memcpy(dest->Row, src->Row, (m + 1) * sizeof(long));
  memcpy(dest->Col, src->Col, (d + 1) * sizeof(long));
}

/*
 * pushQ(lrs_dat *globals,m,d):
 * this routine ensures that Qtail points to a record that
 * may be copied into.
 *
 * It may create a new record, or it may just move the head pointer
 * forward so that know that the old record has been overwritten.
 */

static void pushQ(lrs_dat *global, long m, long d, long m_A) {

  if ((global->Qtail->next) == global->Qhead) {
    /* the Queue is full */
    if (dict_count < dict_limit) {
      /* but we are allowed to create more */
      lrs_dic *p;

      p = new_lrs_dic(m, d, m_A);

      if (p) {

        /* we successfully created another record */

        p->next = global->Qtail->next;
        (global->Qtail->next)->prev = p;
        (global->Qtail->next) = p;
        p->prev = global->Qtail;

        dict_count++;
        global->Qtail = p;

      } else {
        /* virtual memory exhausted. bummer */
        global->Qhead = global->Qhead->next;
        global->Qtail = global->Qtail->next;
      }
    } else {
      /*
       * user defined limit reached. start overwriting the
       * beginning of Q
       */
      global->Qhead = global->Qhead->next;
      global->Qtail = global->Qtail->next;
    }
  }

  else {
    global->Qtail = global->Qtail->next;
  }
}

lrs_dic *lrs_getdic(lrs_dat *Q)
/* create another dictionary for Q without copying any values */
/* derived from lrs_alloc_dic,  used by nash.c                */
{
  lrs_dic *p;

  long m;

  m = Q->m;

  p = new_lrs_dic(m, Q->n - 1, Q->m);
  if (!p)
    return NULL;

  p->next = p;
  p->prev = p;
  Q->Qhead = p;
  Q->Qtail = p;

  return p;
}

#define NULLRETURN(e)                                                          \
  if (!(e))                                                                    \
    return NULL;

static lrs_dic *new_lrs_dic(long m, long d, long m_A) {
  lrs_dic *p;

  NULLRETURN(p = (lrs_dic *)malloc(sizeof(lrs_dic)));

  NULLRETURN(p->B = (long int *)calloc((m + 1), sizeof(long)));
  NULLRETURN(p->Row = (long int *)calloc((m + 1), sizeof(long)));

  NULLRETURN(p->C = (long int *)calloc((d + 1), sizeof(long)));
  NULLRETURN(p->Col = (long int *)calloc((d + 1), sizeof(long)));

#if defined(GMP)
  lrs_alloc_mp(p->det);
  lrs_alloc_mp(p->objnum);
  lrs_alloc_mp(p->objden);
#endif

  p->d_orig = d;
  p->A = lrs_alloc_mp_matrix(m_A, d);

  return p;
}

void lrs_free_dic(lrs_dic *P, lrs_dat *Q) {
  /* do the same steps as for allocation, but backwards */
  /* gmp variables cannot be cleared using free: use lrs_clear_mp* */
  lrs_dic *P1;

  if (Q == NULL) {
    return;
  }

  if (P == NULL) {
    return;
  }
  /* repeat until cache is empty */

  do {
    /* I moved these here because I'm not certain the cached dictionaries
       need to be the same size. Well, it doesn't cost anything to be safe. db
     */

    long d = P->d_orig;
    long m_A = P->m_A;

    lrs_clear_mp_matrix(P->A, m_A, d);

    /* "it is a ghastly error to free something not assigned my malloc" KR167 */
    /* so don't try: free (P->det);                                           */

    lrs_clear_mp(P->det);
    lrs_clear_mp(P->objnum);
    lrs_clear_mp(P->objden);

    free(P->Row);
    free(P->Col);
    free(P->C);
    free(P->B);

    /* go to next record in cache if any */
    P1 = P->next;
    free(P);
    P = P1;

  } while (Q->Qhead != P);

  Q->Qhead = NULL;
  Q->Qtail = NULL;
}

void lrs_free_dat(lrs_dat *Q) {

  int i = 0;

  if (Q == NULL) {
    return;
  }

  /* most of these items were allocated in lrs_alloc_dic */

  lrs_clear_mp_vector(Q->Gcd, Q->m);
  lrs_clear_mp_vector(Q->Lcm, Q->m);

  lrs_clear_mp(Q->Nvolume);

  free(Q->redundcol);
  free(Q->inequality);
  free(Q->linearity);
  free(Q->minratio);
  free(Q->temparray);

  /*2020.8.1 DA: lrs_global_list is not a stack but a list, so have to delete Q
   */

  while (i < lrs_global_count && lrs_global_list[i] != Q)
    i++;

  if (i != lrs_global_count)
    while (i < lrs_global_count) {
      lrs_global_list[i] = lrs_global_list[i + 1];
      i++;
    }

  lrs_global_count--;
  free(Q);
}

static long check_cache(lrs_dic **D_p, lrs_dat *global, long *i_p, long *j_p) {
  /* assign local variables to structures */

  cache_tries++;

  if (global->Qtail == global->Qhead) {
    /* Q has only one element */
    cache_misses++;
    return 0;

  } else {
    global->Qtail = global->Qtail->prev;

    *D_p = global->Qtail;

    *i_p = global->Qtail->i;
    *j_p = global->Qtail->j;

    return 1;
  }
}

lrs_dic *lrs_alloc_dic(lrs_dat *Q)
/* allocate and initialize lrs_dic */
{

  lrs_dic *p;
  long i, j;
  long m, d, m_A;

  m = Q->m;
  d = Q->n - 1;
  m_A = m; /* number of rows in A */

  p = new_lrs_dic(m, d, m_A);
  if (!p)
    return NULL;

  p->next = p;
  p->prev = p;
  Q->Qhead = p;
  Q->Qtail = p;

  dict_count = 1;
  dict_limit = 50;
  cache_tries = 0;
  cache_misses = 0;

  /* Initializations */

  p->d = p->d_orig = d;
  p->m = m;
  p->m_A = m_A;
  p->depth = 0L;
  p->lexflag = TRUE;
  itomp(ONE, p->det);
  itomp(ZERO, p->objnum);
  itomp(ONE, p->objden);

  /*m+d+1 is the number of variables, labelled 0,1,2,...,m+d  */
  /*  initialize array to zero   */
  for (i = 0; i <= m_A; i++)
    for (j = 0; j <= d; j++)
      itomp(ZERO, p->A[i][j]);

  if (Q->nlinearity == ZERO) /* linearity may already be allocated */
    Q->linearity = (long int *)CALLOC((m + d + 1), sizeof(long));

  Q->inequality = (long int *)CALLOC((m + d + 1), sizeof(long));
  Q->redundcol = (long int *)CALLOC((m + d + 1), sizeof(long));
  Q->minratio = (long int *)CALLOC((m + d + 1), sizeof(long));
  Q->temparray = (long int *)CALLOC((unsigned)m + d + 1, sizeof(long));

  Q->inequality[0] = 2L;
  Q->Gcd = lrs_alloc_mp_vector(m);
  Q->Lcm = lrs_alloc_mp_vector(m);

  Q->lastdv = d; /* last decision variable may be decreased */
                 /* if there are redundant columns          */

  for (i = 0; i < m + d + 1; i++) {
    Q->inequality[i] = 0;
  }

  for (i = 0; i <= m; i++) {
    if (i == 0)
      p->B[0] = 0;
    else
      p->B[i] = d + i;
    p->Row[i] = i;
  }

  for (j = 0; j < d; j++) {
    p->C[j] = j + 1;
    p->Col[j] = j + 1;
  }
  p->C[d] = m + d + 1;
  p->Col[d] = 0;
  return p;
} /* end of lrs_alloc_dic */

/* Routines based on lp_solve */

void lrs_set_row(lrs_dic *P, lrs_dat *Q, long row, long num[], long den[],
                 long ineq)
/* convert to lrs_mp then call lrs_set_row */
{
  lrs_mp_vector Num, Den;
  long d;
  long j;

  d = P->d;

  Num = lrs_alloc_mp_vector(d + 1);
  Den = lrs_alloc_mp_vector(d + 1);

  for (j = 0; j <= d; j++) {
    itomp(num[j], Num[j]);
    itomp(den[j], Den[j]);
  }

  lrs_set_row_mp(P, Q, row, Num, Den, ineq);

  lrs_clear_mp_vector(Num, d + 1);
  lrs_clear_mp_vector(Den, d + 1);
}

void lrs_set_row_mp(lrs_dic *P, lrs_dat *Q, long row, lrs_mp_vector num,
                    lrs_mp_vector den, long ineq)
/* set row of dictionary using num and den arrays for rational input */
/* ineq = 1 (GE)   - ordinary row  */
/*      = 0 (EQ)   - linearity     */
{
  lrs_mp Temp, mpone;
  lrs_mp_vector oD; /* denominator for row  */

  long i, j;

  /* assign local variables to structures */

  lrs_mp_matrix A;
  lrs_mp_vector Gcd, Lcm;
  long d;
  lrs_alloc_mp(Temp);
  lrs_alloc_mp(mpone);
  A = P->A;
  d = P->d;
  Gcd = Q->Gcd;
  Lcm = Q->Lcm;

  oD = lrs_alloc_mp_vector(d);
  itomp(ONE, mpone);
  itomp(ONE, oD[0]);

  i = row;
  itomp(ONE, Lcm[i]);         /* Lcm of denominators */
  itomp(ZERO, Gcd[i]);        /* Gcd of numerators */
  for (j = 0; j <= d; j++)
  {
    copy(A[i][j], num[j]);
    copy(oD[j], den[j]);
    if (!one(oD[j]))
      lcm(Lcm[i], oD[j]); /* update lcm of denominators */
    copy(Temp, A[i][j]);
    gcd(Gcd[i], Temp); /* update gcd of numerators   */
  }

  storesign(Gcd[i], POS);
  storesign(Lcm[i], POS);
  if (mp_greater(Gcd[i], mpone) || mp_greater(Lcm[i], mpone))
    for (j = 0; j <= d; j++) {
      exactdivint(A[i][j], Gcd[i], Temp); /*reduce numerators by Gcd  */
      mulint(Lcm[i], Temp, Temp);         /*remove denominators */
      exactdivint(Temp, oD[j], A[i][j]);  /*reduce by former denominator */
    }

  if (ineq == EQ) /* input is linearity */
  {
    Q->linearity[Q->nlinearity] = row;
    Q->nlinearity++;
  }

  lrs_clear_mp_vector(oD, d);
  lrs_clear_mp(Temp);
  lrs_clear_mp(mpone);
} /* end of lrs_set_row_mp */

long dan_selectpivot(lrs_dic *P, lrs_dat *Q, long *r, long *s)
/* select pivot indices using dantzig simplex method             */
/* largest coefficient with lexicographic rule to avoid cycling  */
/* Bohdan Kaluzny's handiwork                                    */
/* returns TRUE if pivot found else FALSE                        */
/* pivot variables are B[*r] C[*s] in locations Row[*r] Col[*s]  */
{
  long j, k, col;
  lrs_mp coeff;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *Col = P->Col;
  long d = P->d;

  /*  printf("\n*dantzig"); */
  lrs_alloc_mp(coeff);
  *r = 0;
  *s = d;
  j = 0;
  k = 0;

  itomp(0, coeff);
  /*find positive cost coef */
  while (k < d) {
    if (mp_greater(A[0][Col[k]], coeff)) {
      j = k;
      copy(coeff, A[0][Col[j]]);
    }
    k++;
  }

  if (positive(coeff)) /* pivot column found! */
  {
    *s = j;
    col = Col[j];

    /*find min index ratio */
    *r = lrs_ratio(P, Q, col);
    if (*r != 0) {
      lrs_clear_mp(coeff);
      return (TRUE); /* pivot found */
    }
  }
  lrs_clear_mp(coeff);
  return (FALSE);
} /* end of dan_selectpivot        */

long ran_selectpivot(lrs_dic *P, lrs_dat *Q, long *r, long *s)
/* select pivot indices using random edge rule                   */
/* largest coefficient with lexicographic rule to avoid cycling  */
/* pivot variables are B[*r] C[*s] in locations Row[*r] Col[*s]  */
{
  long i, j, k, col, t;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *Col = P->Col;
  long d = P->d;
  long *perm;

  perm = (long *)calloc((d + 1), sizeof(long));
  *r = 0;
  *s = d;
  k = 0;
  /*  printf("\n*random edge"); */

  /* generate random permutation of 0..d-1 */
  for (i = 0; i < d; i++)
    perm[i] = i;

  for (i = 0; i < d; i++) {
    j = rand() % (d - i) + i;
    t = perm[j];
    perm[j] = perm[i];
    perm[i] = t; // Swap i and j
  }

  /*find first positive cost coef according to perm */
  while (k < d && !positive(A[0][Col[perm[k]]]))
    k++;

  if (k < d) /* pivot column found! */
  {
    j = perm[k];
    *s = j;
    col = Col[j];

    /*find min index ratio */
    *r = lrs_ratio(P, Q, col);
    if (*r != 0) {
      free(perm);
      return (TRUE); /* pivot found */
    }
  }
  free(perm);
  return (FALSE);
} /* end of ran_selectpivot        */

void lrs_overflow(int parm) {
  printf("lrs_overflow called");
  lrs_exit(parm); 
  }

void lrs_exit(int i) {
  // fflush(stdout);
  // exit(i);
  printf("lrsexit called (maybe overflow)");
  throw std::exception();
}
