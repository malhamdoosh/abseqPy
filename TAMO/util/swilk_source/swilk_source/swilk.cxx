/* SWILK -- Compute Shapiro-Wilk test metric and p-value */

/* Author: David Benjamin Gordon (hand-ported from F77) */

/* C-code Copyright 2005, Whitehead Institute for Biomedical Research.  All rights reserved. */

/* Original F77 Code from Statlib (Applied Statistics) by Patrick Royston                    */ 

/* Copyright (2005) Whitehead Institute for Biomedical Research */
/* All Rights Reserved */

#include <stdio.h>
#include <math.h>


void  test();
void  swilkf(int &init, float *x, int &n, int &n1, int &n2, float *a, float &w, float &pw, int &ifault);
float poly(float *C, int NORD, float X);
float ppnd(float p);
float sndev(float PROB);
float alnorm(float x, int upper);
void  stat8(float Z, float &Q, float &P);
float min(float a, float b);
float pow(float a, float b);
float sign(int a, float b);


/*
c     AS R94 -- calculates Shapiro-Wilk normality test and P-value
c     for sample sizes 3 <= n <= 5000 .  Handles censored or uncensored data.
c     Corrects AS 181, which was found to be inaccurate for n > 50.
c     by Patrick Royston, Applied Statistics vol. 44 no. 4 (1995)
c
c     Example driver routine for SWILK.  Includes data from original
c     article.  Published results for ncens=0:  w=.83467, pw=.000914
c
*/

void test() {
  int   n = 25;
  float a[25];
  float x[25] = {.139,.157,.175,.256,.344,.413,.503,.577,.614,.655,
		 .954,1.392,1.557,1.648,1.690,1.994,2.174,2.206,3.245,3.510,
		 3.571,4.354,4.980,6.084,8.351};
  /*
      parameter (n=25)
      real x(n),a(n)
      logical init
      number censored
      ncens=0
      n1=n-ncens
      n2=n/2
      call swilk(init, x, n, n1, n2, a, w, pw, ifault)
      write(*,101)w,pw,ifault
101   format(' w = ',f9.6,' pw = ',f10.6,' ifault = ',i1)
      stop
      end
  */
  for (int i=0;i<25;i++) a[i]=0.0;
  int init = 0;
  int ncens = 0;
  int n1 = n - ncens;
  int n2 = n1/2;
  float w,pw;
  int ifault;

  w = 0;
  pw = 0;


  swilkf(     init,        x,      n,      n1,      n2,        a,        w,        pw,      ifault);
  

  printf(" w = %9.6f    pw = %10.6f   ifault = %d\n",w,pw,ifault);
}

/*
      subroutine swilk(init, x, n, n1, n2, a, w, pw, ifault)
c
c     AS R94  Calculates the Shapiro-Wilk W test and its significance level.
c     by Patrick Royston, Applied Statistics vol. 44 no. 4 (1995)
c
c      Arguments for swilk:
c      Inputs:
c      ~~~~~~~
c      init:   flag, .false. to calc coefficients a() for W test.
c      x(n1):  input data, must be sorted in ascending order
c      n:      (uncensored) sample size
c      n1:     number of uncensored values (n1=n for a complete sample)
c      n2:     n/2
c      a(n2):  (if init is .true.), coefficients for W test
c              (if init is .false.), undefined
c      w:      zero or positive to calculate W from x(n1),
c              negative to calculate P-value of -w.
c      Outputs:
c      ~~~~~~~~
c      init:   .true. unless an error has occurred (ifault.gt.0)
c      a(n2):  coefficients for W test
c      w:      W test statistic
c      pw:     P-value for W
c      ifault: fault indicator, all non-zero codes fatal except ifault.eq.2:
c            0: ok
c            1: n1 < 3
c            2: n > 5000 (non-fatal error)
c            3: n2 < n/2, so insufficient storage for a()
c            4: n1 > n or ((n1 < n) & (n < 20))
c            5: proportion censored, (n-n1)/n, > 0.8
c            6: effectively zero range (assuming data sorted)
c            7: data not in increasing sort order
c
      integer n, n1, n2, ifault
      real x(*), a(*), pw, w
*/
void swilkf(int &init, float *x, int &n, int &n1, int &n2, float *a, float &w, float &pw, int &ifault) {
  /*
    real c1(6), c2(6), c3(4), c4(4), c5(4), c6(3), c7(2)
      real c8(2), c9(2), g(2)
      real z90, z95, z99, zm, zss, bf1, xx90, xx95, zero, one, two
      real three, sqrth, qtr, th, small, pi6, stqr
c
      logical init, upper
  */
  /*
      data c1/0.0e0, 0.221157e0, -0.147981e0, -0.207119e1,
     +       0.4434685e1, -0.2706056e1/
      data c2/0.0e0, 0.42981e-1, -0.293762e0, -0.1752461e1,
     +       0.5682633e1, -0.3582633e1/
      data c3/0.5440e0, -0.39978e0, 0.25054e-1, -0.6714e-3/
      data c4/0.13822e1, -0.77857e0, 0.62767e-1, -0.20322e-2/
      data c5/-0.15861e1, -0.31082e0, -0.83751e-1, 0.38915e-2/
      data c6/-0.4803e0, -0.82676e-1, 0.30302e-2/
      data c7/0.164e0, 0.533e0/
      data c8/0.1736e0, 0.315e0/
      data c9/0.256e0, -0.635e-2/
      data g/-0.2273e1, 0.459e0/
      data z90, z95, z99/0.12816e1, 0.16449e1, 0.23263e1/
      data zm, zss/0.17509e1, 0.56268e0/
      data bf1/0.8378e0/, xx90, xx95/0.556e0, 0.622e0/
      data zero/0.0e0/, one/1.0e0/, two/2.0e0/, three/3.0e0/
      data sqrth/0.70711e0/, qtr/0.25e0/, th/0.375e0/, small/1e-19/
      data pi6/0.1909859e1/, stqr/0.1047198e1/, upper/.true./
  */
  float c1[] = {0.0e0, 0.221157e0, -0.147981e0, -0.207119e1, 0.4434685e1, -0.2706056e1};
  float c2[] = {0.0e0, 0.42981e-1, -0.293762e0, -0.1752461e1,  0.5682633e1, -0.3582633e1};
  float c3[] = {0.5440e0, -0.39978e0, 0.25054e-1, -0.6714e-3};
  float c4[] = {0.13822e1, -0.77857e0, 0.62767e-1, -0.20322e-2};
  float c5[] = {-0.15861e1, -0.31082e0, -0.83751e-1, 0.38915e-2};
  float c6[] = {-0.4803e0, -0.82676e-1, 0.30302e-2};

  float c7[] = {0.164e0, 0.533e0};
  float c8[] = {0.1736e0, 0.315e0};
  float c9[] = {0.256e0, -0.635e-2};
  float g[] = {-0.2273e1, 0.459e0};
  float z90   = 0.12816e1;
  float z95   = 0.16449e1;
  float z99   = 0.23263e1;
  float zm    = 0.17509e1;
  float zss   = 0.56268e0;
  float bf1   = 0.8378e0;
  float xx90  = 0.556e0;
  float xx95  = 0.622e0;
  float zero = 0.0;
  float one = 1.0;
  float two = 2.0;
  float three = 3.0;
  float sqrth = 0.70711e0;
  float qtr   = 0.25;
  float th    = 0.375;
  float small = 1e-19;
  float pi6   = 0.1909859e1;
  float stqr = 0.1047198e1;
  int   upper = 1;

  /*
      real summ2, ssumm2, fac, rsn, an, an25, a1, a2, delta, range
      real sa, sx, ssx, ssa, sax, asa, xsx, ssassx, w1, y, xx, xi
      real gamma, m, s, ld, bf, z90f, z95f, z99f, zfm, zsd, zbar
      real ppnd, alnorm, poly
      integer ncens, nn2, i, i1, j
  */
  float summ2, ssumm2, fac, rsn, an, an25, a1, a2, delta, range;
  float sa, sx, ssx, ssa, sax, asa, xsx, ssassx, w1, y, xx, xi;
  float gamma, m, s, ld, bf, z90f, z95f, z99f, zfm, zsd, zbar;
  int   ncens, nn2, i, i1, j;

  /*
      pw=one
      if (w.ge.zero) w=one
      an=n
      ifault=3
      nn2=n/2
      if (n2.lt.nn2) return
      ifault=1
      if (n.lt.3) return
  */
  pw=one;
  if (w >= zero) w=one;
  an=n;
  ifault=3;
  nn2=n/2;
  if (n2 < nn2) return;
  ifault=1;
  if (n < 3) return;
	/*
c
c      If init is false, calculate coefficients for the test.
c
      if (.not.init) then
            if (n.eq.3) then
                  a(1)=sqrth
            else
                  an25=an+qtr
                  summ2=zero
                  do 30 i=1,n2
                        a(i)=ppnd((i-th)/an25)
                        summ2=summ2+a(i)**2
30                continue
                  summ2=summ2*two
                  ssumm2=sqrt(summ2)
                  rsn=one/sqrt(an)
                  a1=poly(c1,6,rsn)-a(1)/ssumm2
	*/
  if (!init) {
    if (n == 3) {
      a[0] = sqrth;
    } else {
      an25 = an + qtr;
      summ2 = zero;
      for (i=0; i< n2; i++) {
	a[i] = ppnd((i+1-th)/an25);
	summ2 += a[i]*a[i];
      }
      summ2  = summ2*two;
      ssumm2 = sqrt(summ2);
      rsn    = one / sqrt(an);
      a1     = poly(c1,6,rsn) - (a[0]/ssumm2);
    /*
c
c      Normalize coefficients
c
                  if (n.gt.5) then
                        i1=3
                        a2=-a(2)/ssumm2+poly(c2,6,rsn)
                        fac=sqrt((summ2-two*a(1)**2-two*
     +                         a(2)**2)/(one-two*a1**2-two*a2**2))
                        a(1)=a1
                        a(2)=a2
                  else
                        i1=2
                        fac=sqrt((summ2-two*a(1)**2)/
     +                         (one-two*a1**2))
                        a(1)=a1
                  end if
                  do 40 i=i1,nn2
40                      a(i)=-a(i)/fac
            end if
            init=.true.
      end if
      if (n1.lt.3) return
      ncens=n-n1
      ifault=4
      if (ncens.lt.0 .or. (ncens.gt.0 .and. n.lt.20)) return
      ifault=5
      delta=float(ncens)/an
      if (delta.gt.0.8) return
    */
      if ( n > 5 ) {
	i1  = 3;
	a2  = -a[1]/ssumm2+poly(c2,6,rsn);
	fac = sqrt((summ2 - two*a[0]*a[0]-two*a[1]*a[1]) / 
		   (one - two*a1*a1 - two*a2*a2));
	a[0] = a1;
	a[1] = a2;
      } else {
	i1  = 2;
	fac = sqrt((summ2-two*a[0]*a[0]) / 
		   (one - two*a1*a1));
	a[0]= a1;
      }
      for (i=i1-1; i<nn2; i++) {
	a[i] = -a[i]/fac;
      } 
    }
    init = 1;
  }
  if (n1 < 3) return;
  ncens = n - n1;
  ifault = 4;
  if ((ncens < 0) || (ncens > 0 and n < 20)) return;
  ifault = 5;
  delta = float(ncens)/an;
  if (delta > 0.8) return;
/*
c
c      If w input as negative, calculate significance level of -w.
c
      if (w.lt.zero) then
            w1=one+w
            ifault=0
            goto 70
      end if
*/
  if (w < 0) {
    w1 = one + w;
    ifault = 0;
    goto __seventy;
  }
/*
c
c      Check for "zero" range
c
      ifault=6
      range=x(n1)-x(1)
      if (range.lt.small) return
*/
  ifault = 6;
  range = x[n1-1] - x[0];
  if (range < small) return;
/*
c
c      Check for correct sort order on range-scaled x
c
      ifault=7
      xx=x(1)/range
      sx=xx
      sa=-a(1)
      j=n-1
      do 50 i=2,n1
            xi=x(i)/range
            if (xx-xi.gt.small) return
            sx=sx+xi
            if (i.ne.j) sa=sa+sign(1,i-j)*a(min(i,j))
            xx=xi
            j=j-1
50    continue
      ifault=0
      if (n.gt.5000) ifault=2
*/
  ifault = 7;
  xx = x[0]/range;
  sx = xx;
  sa = -a[0];
  j  = n-1-1   ; // NEED EXTRA -1?
  for (i=1; i<n1; i++) {
    xi = x[i]/range;
    if (xx-xi > small) return;
    sx = sx + xi;
    if (i != j) sa = sa + sign(1,i-j) * a[int(min(i,j))];
    xx = xi;
    j=j-1;
  }
  ifault = 0;
  if (n > 5000) ifault = 2;
  
/*
c
c      Calculate W statistic as squared correlation
c      between data and coefficients.
c
      sa=sa/n1
      sx=sx/n1
      ssa=zero
      ssx=zero
      sax=zero
      j=n
      do 60 i=1,n1
            if (i.ne.j) then
                  asa=sign(1,i-j)*a(min(i,j))-sa
            else
                  asa=-sa
            end if
            xsx=x(i)/range-sx
            ssa=ssa+asa*asa
            ssx=ssx+xsx*xsx
            sax=sax+asa*xsx
            j=j-1
60    continue
c
c      w1 equals 1-w calculated to avoid excessive rounding error
c      for w very near 1 -- a potential problem in very large samples
c
      ssassx=sqrt(ssa*ssx)
      w1=(ssassx-sax)*(ssassx+sax)/(ssa*ssx)
*/
  sa=sa/n1;
  sx=sx/n1;
  ssa=zero;
  ssx=zero;
  sax=zero;
  j=n-1; // SHOULD THIS BE n-1?
  for (i=0; i< n1; i++) {
    if (i!=j) {
      asa = sign(1,i-j) * a[int(min(i,j))] - sa;
    } else {
      asa = -sa;
    }
    xsx=x[i]/range-sx;
    ssa=ssa+asa*asa;
    ssx=ssx+xsx*xsx;
    sax=sax+asa*xsx;
    j=j-1;
  }
  ssassx=sqrt(ssa*ssx);
  w1=(ssassx-sax)*(ssassx+sax)/(ssa*ssx);

/*
70    w=one-w1
*/
 __seventy:
  w = one-w1;
/*
c
c      Calculate significance level for W (exact for n=3).
c
      if (n.eq.3) then
            pw=pi6*(asin(sqrt(w))-stqr)
            return
      end if
      y=log(w1)
      xx=log(an)
      m=zero
      s=one
      if (n.le.11) then
            gamma=poly(g,2,an)
            if (y.ge.gamma) then
                  pw=small
                  return
            end if
            y=-log(gamma-y)
            m=poly(c3,4,an)
            s=exp(poly(c4,4,an))
      else
            m=poly(c5,4,xx)
            s=exp(poly(c6,3,xx))
      end if
*/
  if (n == 3) { 
    pw = pi6 * asin(sqrt(w))-stqr;
    return;
  } 
  y = log(w1);
  xx = log(an);
  m=zero;
  s=one;
  if (n <= 11) {
    gamma=poly(g,2,an);
    if (y >= gamma) { 
      pw=small;
      return;
    }
    y=-log(gamma-y);
    m=poly(c3,4,an);
    s=exp(poly(c4,4,an));
  } else {
    m=poly(c5,4,xx);
    s=exp(poly(c6,3,xx));
  }
/*
c
c      Censoring by proportion ncens/n.
c      Calculate mean and sd of normal equivalent deviate of W.
c
      if (ncens.gt.0) then
            ld=-log(delta)
            bf=one+xx*bf1
            z90f=z90+bf*poly(c7,2,xx90**xx)**ld
            z95f=z95+bf*poly(c8,2,xx95**xx)**ld
            z99f=z99+bf*poly(c9,2,xx)**ld
c
c      Regress z90f,...,z99f on normal deviates z90,...,z99 to get
c      pseudo-mean and      pseudo-sd of z as the slope and intercept.
c
            zfm=(z90f+z95f+z99f)/three
            zsd=(z90*(z90f-zfm)+z95*(z95f-zfm)+z99*(z99f-zfm))/zss
            zbar=zfm-zsd*zm
            m=m+zbar*s
            s=s*zsd
      end if
      pw=alnorm((y-m)/s,upper)
      return
      end
*/
  if (ncens > 0 ) {
    ld = -log(delta);
    bf=one+xx*bf1;
    z90f=z90+bf*pow(poly(c7,2,pow(xx90,xx)),ld);
    z95f=z95+bf*pow(poly(c8,2,pow(xx95,xx)),ld);
    z99f=z99+bf*pow(poly(c9,2,xx),ld);
    
    zfm=(z90f+z95f+z99f)/three;
    zsd=(z90*(z90f-zfm)+z95*(z95f-zfm)+z99*(z99f-zfm))/zss;
    zbar=zfm-zsd*zm;
    m=m+zbar*s;
    s=s*zsd;
  }

  pw = alnorm((y-m)/s,upper);
}

/*
c
c-----------------------------------------------------------------------
c      Auxiliary routines
c-----------------------------------------------------------------------
      FUNCTION POLY(C,NORD,X)
C
C          CALCULATES THE ALGEBRAIC POLYNOMIAL OF ORDER NORD-1 WITH
C          ARRAY OF COEFFICIENTS C. ZERO ORDER COEFFICIENT IS C(1).
C
*/
float poly(float *C, int NORD, float X) {
  /*
      DIMENSION C(NORD)
      POLY=C(1)
      IF(NORD.EQ.1)RETURN
      P=X*C(NORD)
      IF(NORD.EQ.2)GOTO20
      N2=NORD-2
      J=N2+1
      DO 10 I=1,N2
      P=(P+C(J))*X
      J=J-1
  10  CONTINUE
  20  POLY=POLY+P
      RETURN
      END
  */
  int i;
  float P;
  float ans;
  int   J, N2;
  

  ans = C[0];
  if (NORD == 1) return ans;
  P = X * C[NORD-1];
  if (NORD == 2) goto __twenty;
  N2 = NORD-2;
  J  = N2;  /* No +1 because of C vs. Foratran nonsense */
  for (i=0; i<N2; i++) {
    P = (P + C[J]) * X;
    J = J-1;
  }
 __twenty:
  ans = ans + P;
  return ans;
}

/*
c
      REAL FUNCTION PPND(P)
      PPND=SNDEV(P)
      RETURN
      END
*/
float ppnd(float p) {
  return sndev(p);
}

/*
c
      real function SNDEV(PROB)
c
c     Modified 30-Jul-91 so if PROB is zero, returns BIG as 99.9999 not 1E+38.
c     Some logic modified to F77 standard.
c
*/
float sndev(float PROB) {
  /*
      DOUBLE PRECISION A0,A1,A2,A3,B1,B2,B3,B4,C0,C1,C2,C3,D1,D2,
     +Q,R,ONE,HALF,SPLIT,ZERO
C
C      BASED ON ALGORITHM AS 111
C
      data big/0.999999e+02/
      DATA ONE/1.0D0/,HALF/0.5D0/,SPLIT/0.42D0/,ZERO/0.0D0/
      DATA A0 /  2.50662823884D0/
      DATA A1 /-18.61500062529D0/
      DATA A2 / 41.39119773534D0/
      DATA A3 /-25.44106049637D0/
      DATA B1 / -8.47351093090D0/
      DATA B2 / 23.08336743743D0/
      DATA B3 /-21.06224101826D0/
      DATA B4 /  3.13082909833D0/
      DATA C0 / -2.78718931138D0/
      DATA C1 / -2.29796479134D0/
      DATA C2 /  4.85014127135D0/
      DATA C3 /  2.32121276858D0/
      DATA D1 /  3.54388924762D0/
      DATA D2 /  1.63706781897D0/
  */
  double Q, R;
  double big = 0.999999e+02;
  double ONE = 1.0e0;
  double HALF = 0.5e0;
  double SPLIT = 0.42e0;
  double ZERO = 0.0e0;
  double A0  =   2.50662823884e0;
  double A1  = -18.61500062529e0;
  double A2  =  41.39119773534e0;
  double A3  = -25.44106049637e0;
  double B1  =  -8.47351093090e0;
  double B2  =  23.08336743743e0;
  double B3  = -21.06224101826e0;
  double B4  =   3.13082909833e0;
  double C0  =  -2.78718931138e0;
  double C1  =  -2.29796479134e0;
  double C2  =   4.85014127135e0;
  double C3  =   2.32121276858e0;
  double D1  =   3.54388924762e0;
  double D2  =   1.63706781897e0;
  /*
      Q=PROB-HALF
      IF(DABS(Q).GT.SPLIT)GOTO 1
      R=Q*Q
      SNDEV=Q*(((A3*R+A2)*R+A1)*R+A0)/
     +  ((((B4*R+B3)*R+B2)*R+B1)*R+ONE)
      RETURN
    1 R=PROB
      IF(Q.GT.ZERO)R=ONE-PROB
      IF(R.GT.ZERO) then
            R=DSQRT(-DLOG(R))
            SNDEV=(((C3*R+C2)*R+C1)*R+C0)/((D2*R+D1)*R+ONE)
      else
            SNDEV=big
      end if
      IF(Q.LT.ZERO)SNDEV=-SNDEV
      RETURN
      END
  */
  float SNDEV; 

  Q=PROB-HALF;
  if (fabs(Q) > SPLIT) goto __one;
  R=Q*Q;
  SNDEV=Q*(((A3*R+A2)*R+A1)*R+A0)/((((B4*R+B3)*R+B2)*R+B1)*R+ONE);
  return SNDEV;
  __one:
  R=PROB;
  if (Q > ZERO) R = ONE-PROB;
  if (R > ZERO) {
    R=sqrt(-log(R));
    SNDEV=(((C3*R+C2)*R+C1)*R+C0)/((D2*R+D1)*R+ONE);
  } else {
    SNDEV=big;
  }
  if (Q < ZERO) SNDEV = -SNDEV;
  return SNDEV;
}
  /*
      FUNCTION ALNORM(X,UPPER)
      LOGICAL UPPER
      CALL STAT8(X,Q,P)
      ALNORM=Q
      IF(UPPER)ALNORM=P
      RETURN
      END
  */
float alnorm(float x, int upper) {
  float Q, P;
  stat8(x, Q, P);
  if (upper) {
    return P;
  } else {
    return Q;
  }
}


/*
c
      SUBROUTINE STAT8(Z,Q,P)
*/
void stat8(float Z, float &Q, float &P) {
/*
C   REAL VERSION
C  Q AND P ARE THE LOWER AND UPPER TAIL AREAS OF A UNIT NORMAL DISTRI-
C  -BUTION AT A DEVIATION X.  THE RESULTS ARE GOOD TO 9S FOR ANY X.  SEE
C  ADAMS,A.G., COMP. J., 12, 197-198, 1969.
C
      REAL Z,Q,P
      DOUBLE PRECISION X,Y,T,HALF,ONE,SPLIT,A0,A1,A2,A3,A4,A5,A6,
     +A7,B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,ZERO,SMALL
      DATA HALF/0.5D0/,ONE/1.0D0/,SPLIT/1.28D0/,SMALL/1.0D-38/
      DATA ZERO/0.0D0/
      DATA A0/3.98942280444D-1/,A1/-3.99903438504D-1/
      DATA A2/5.75885480458D0/,A3/-2.98213557808D+1/
      DATA A4/2.62433121679D0/,A5/4.86959930692D+1/
      DATA A6/5.92885724438D0/,A7/3.98942280385D-1/
      DATA B0/-3.8052D-8/,B1/1.00000615302D0/
      DATA B2/3.98064794D-4/,B3/1.98615381364D0/
      DATA B4/-1.51679116635D-1/,B5/5.29330324926D0/
      DATA B6/4.8385912808D0/,B7/-1.51508972451D+1/
      DATA B8/7.42380924027D-1/,B9/3.0789933034D+1/
      DATA B10/3.99019417011D0/
*/
  double X, Y, T;
  double L;
  double HALF = 0.5e0;
  double ONE = 1.0e0;
  double SPLIT = 1.28e0;
  double SMALL = 1.0e-38;
  double ZERO = 0.0e0;
  double A0 = 3.98942280444e-1;
  double A1 = -3.99903438504e-1;
  double A2 = 5.75885480458e0;
  double A3 = -2.98213557808e+1;
  double A4 = 2.62433121679e0;
  double A5 = 4.86959930692e+1;
  double A6 = 5.92885724438e0;
  double A7 = 3.98942280385e-1;
  double B0 = -3.8052e-8;
  double B1 = 1.00000615302e0;
  double B2 = 3.98064794e-4;
  double B3 = 1.98615381364e0;
  double B4 = -1.51679116635e-1;
  double B5 = 5.29330324926e0;
  double B6 = 4.8385912808e0;
  double B7 = -1.51508972451e+1;
  double B8 = 7.42380924027e-1;
  double B9 = 3.0789933034e+1;
  double B10 = 3.99019417011e0;


/*
      IF (ABS(Z).GE.13.0) GO TO 18
      X=Z
      L=0
      IF(X)10,11,12
10    L=1
      X=-X
  12  Y=HALF*X*X
      IF(X.GE.SPLIT)GOTO13
      T=HALF-X*(A0+A1*Y/(Y+A2+A3/(Y+A4+A5/(Y+A6))))
      GOTO14
  13  Y=A7*DEXP(-Y)
      IF(DABS(Y/X).GT.SMALL)GOTO15
      T=ZERO
      GOTO14
  15  T=Y/(X+B0+B1/(X+B2+B3/(X+B4+B5/(X+B6+B7/(X+B8+B9/(X+B10))))))
14    IF(L.EQ.0)GOTO17
      Q=T
      P=ONE-T
      RETURN
17    P=T
      Q=ONE-T
      RETURN
   11 P=HALF
      Q=HALF
      RETURN
   18 CONTINUE
      IF (Z) 19,11,20
   19 P = 1.0
      Q = 0.0
      RETURN
   20 P = 0.0
      Q = 1.0
      RETURN
      END
*/
  
  if (fabs(Z) > 13.0) goto __L18;
  X = Z;
  L = 0;
  if      (X < 0.0) goto __L10;
  else if (X ==0.0) goto __L11;
  else           goto __L12;
 __L10:
  L=1;
  X=-X;
 __L12:
  Y=HALF*X*X;
  if (X >= SPLIT) goto __L13;
  T=HALF-X*(A0+A1*Y/(Y+A2+A3/(Y+A4+A5/(Y+A6))));
  goto __L14;
 __L13:
  Y=A7*exp(-Y);
  if (fabs(Y/X) > SMALL) goto __L15;
  T=ZERO;
  goto __L14;
 __L15:
  T=Y/(X+B0+B1/(X+B2+B3/(X+B4+B5/(X+B6+B7/(X+B8+B9/(X+B10))))));
 __L14:
  if (L == 0.0) goto __L17;
  Q=T;
  P=ONE-T;
  return;
 __L17:
  P=T;
  Q=ONE-T;
  return;
 __L11:
  P=HALF;
  Q=HALF;
  return;
 __L18:
  if      (Z < 0.0) goto __L19;
  else if (Z ==0.0) goto __L11;
  else              goto __L20;
 __L19:
  P = 1.0;
  Q = 0.0;
  return;
 __L20:
  P = 0.0;
  Q = 1.0;
  return;
}


float min(float a, float b) {
  if (a < b) return a;
  else       return b;
}

float pow(float a, float b) {
  /* Return a**b */
  return exp(b * log(a));
}

float sign(int a, float b) {
  /* A function that assigns the sign of the second argument to the
     absolute value of the first.
  */
  if (b >= 0) return a;
  else        return -a;
}
