#include "MDsupport.h"
#include <stdio.h>
#include <math.h>

// Copyright (2005) Whitehead Institute for Biomedical Research
// All Rights Reserved
// Author: David Benjamin Gordon



int *seq2int(const char *seq, int len) {
  int i;
  if (len==0) for (len=0; seq[len]; len++);
  int *ans = new int[len];
  for (i=0; i<len; i++) {
    switch (seq[i]) {
      case 'A': ans[i] = 0; break;
      case 'C': ans[i] = 1; break;
      case 'G': ans[i] = 2; break;
      case 'T': ans[i] = 3; break;
      case 'N': ans[i] = 4; break;
      default:
	/* printf("# Strange character %c at sequence position %d. Considering as \"N\"\n",seq[i],i); */
	ans[i] = 4;
    }
  }
  return(ans);
}

double *list2double(dlist dptr) {
  return(dptr);
}

void printdouble(double *D, int len) {
  for (int i=0; i<len; i++) printf("%4d: %f\n",i,D[i]);
}


void print_seq(int *iseq,int len) {
  int i;
  for (i=0;i<len;i++) printf("%d",iseq[i]);
  printf("\n");
}

float log2_sum (float logx,float logy) {
  float loga, logb, ans;

  if (logx > logy) { 
    loga = logx;  logb = logy;
  } else {
    logb = logx;  loga = logy;
  }

  if  (logb < -1e100) {
    ans = loga;
  } else {
    ans = loga + log(1 + pow(2,(logb - loga))) * 1.44269504088896 ; /* *1.44 = /log(2) */
  }
  return (ans);
}

PyObject* tupleize_lists(float *probs, float *Zs, int length) {
  int i;
  PyObject *probList = PyList_New(length);
  PyObject *ZList    = PyList_New(length);
  for (i=0; i<length; i++) {
    PyList_SetItem(probList, i, PyFloat_FromDouble(probs[i]));
    PyList_SetItem(ZList,    i, PyFloat_FromDouble(   Zs[i]));
  }
  
  PyObject *PyResult  = PyTuple_New(2);
  PyTuple_SetItem(PyResult, 0, probList);
  PyTuple_SetItem(PyResult, 1,    ZList);
  
  return(PyResult);
  
}


void SeqMat::_init_(int w) {
  if (w < 0 || w>this->_W) {
    printf("Width out of bounds (%d)\n",w);
    return;
  }
  this->width = w;
  for (int i=0; i<w; i++) {
    for (int j=0; j<this->_H; j++) {
      this->M[i][j] = 0L;
    }
    this->mask[i] = 1.0;
  }
  this->beta     = 0.001;
  this->gamma    = 0.5;
  this->gamma_wt = 0.5;
  this->deltamin = 0.001;
}

void SeqMat::compute_ambig() {
}

flist SeqMat::score_probe(char *seq) {
  int i,j;
  int len;
  int *iseq  = seq2int(seq);
  int w;

  w = this->width;
  for (len=0; seq[len]; len++);
  flist ans = new float[len+1];
  for (i=0; i<=len-w; i++) {
    ans[i] = 0;
    /* printf("\n"); */
    for (j=0;j<w;j++) {
      /* printf("%d(%d):%f\t",j,iseq[i+j],this->M[j][iseq[i+j]]); */
      ans[i] += this->M[j][iseq[i+j]];
    }
  }
  ans[i] = -1000;
  return(ans);
}

double SeqMat::scanbest(const char *seq) {
  double best = -100000;
  double tf,tr;
  int i,offset;
  int w = this->width;
  int *iseq;
  int len;

  for (len=0; seq[len]; len++);
  iseq = seq2int(seq,len);

  int comp[] = {3, 2, 1, 0, 4};

  for (offset = 0; offset+w <= len; offset++) {
    tf = 0;
    tr = 0;
    for (i=0; i<w; i++) {
      tf += this->M[i][iseq[offset+i]];
      tr += this->M[w-i-1][comp[iseq[offset+i]]];
      /* tr += this->M[w-i-1][3-iseq[offset+i]]; */
    }
    if (tf>best) best=tf;
    if (tr>best) best=tr;
  } 
  delete[] iseq;
  return(best);    
}

double SeqMat::sumscoresabove(char *seq, float cutoff) {
  float sum  = 0;
  float tf,tr,score;
  int i,offset;
  int w = this->width;
  int *iseq;
  int len;

  for (len=0; seq[len]; len++);
  iseq = seq2int(seq,len);

  int comp[] = {3, 2, 1, 0, 4};

  for (offset = 0; offset+w <= len; offset++) {
    tf = 0;
    tr = 0;
    for (i=0; i<w; i++) {
      tf += this->M[i][iseq[offset+i]];
      tr += this->M[w-i-1][comp[iseq[offset+i]]];
      /* tr += this->M[w-i-1][3-iseq[offset+i]]; */
    }
    if (tf>tr) score = tf;
    else       score = tr;
    if (score >= cutoff) {
      sum += score;
      offset += w-1;
    }
  } 
  delete[] iseq;
  return(sum);    
}

std::vector<int> SeqMat::matchstarts(char *seq, float cutoff) {
  double tf,tr;
  int i,offset;
  int w = this->width;
  int *iseq;
  int len;

  std::vector<int> ans;
  ans.clear();

  for (len=0; seq[len]; len++);
  iseq = seq2int(seq,len);

  int comp[] = {3, 2, 1, 0, 4};

  for (offset = 0; offset+w <= len; offset++) {
    tf = 0;
    tr = 0;
    for (i=0; i<w; i++) {
      tf += this->M[i][iseq[offset+i]];
      tr += this->M[w-i-1][comp[iseq[offset+i]]];
      /* tr += this->M[w-i-1][3-iseq[offset+i]]; */
    }
    if (tf >= cutoff || tr >= cutoff) {
      ans.push_back(offset);
    }
  } 
  delete[] iseq;
  return(ans);
}


double SeqMat::score(int *idxs) {
  double t = 0;
  int i;
  for (i=0; i<this->width; i++) {
    t += this->M[i][idxs[i]];
  }
  return(t);
}

void SeqMat::set(int i, int j, double value) {
  if (i+1> this->width) {
    printf("%d out of bounds (>%d)\n",i,width);
    return;
  }
  if (j+1> this->_H) {
    printf("%d out of bounds (>%d)\n",j,_H);
    return;
  }
  this->M[i][j] = value;
}

double SeqMat::get(int i, int j) {
  if (i+1> this->width) {
    printf("%d out of bounds (>%d)\n",i,width);
    return 0;
  }
  if (j+1> this->_H) {
    printf("%d out of bounds (>%d)\n",j,_H);
    return 0;
  }
  return(this->M[i][j]);
}

double SeqMat::get_c(int i, int j) {
  if (i+1> this->width) {
    printf("%d out of bounds (>%d)\n",i,width);
    return 0;
  }
  if (j+1> this->_H) {
    printf("%d out of bounds (>%d)\n",j,_H);
    return 0;
  }
  return(this->counts[i][j]);
}

void SeqMat::setBg(double A, double C, double G, double T) {
  this->bg[0] = A; this->bg[1] = C; this->bg[2] = G; this->bg[3] = T;
}

void SeqMat::setmask(int pos, float f) {
  if (pos >= this->width || f < 0 || f > 1) {
    printf("Error setting mask pos %d with %f (width %d)\n",
	   pos,f,this->width);
  } else {
    this->mask[pos] = f;
  }
}


void SeqMat::zoops_e(int *iseq, int iseqlen,
		   double *wmerbgs, 
		   float probebg, float gamma,
		   float *probs, float *Zs) {  /* For output */
  int i,j;
  int w           =  this->width;
  int nwmers      =  iseqlen - w + 1;
  float sumScores = -1e101;
  
  /* Compute individual probabilities and total probability */
  /* Primary strand */
  for (i=0; i<nwmers; i++) {
    probs[i] = wmerbgs[i];
    probs[nwmers+i] = wmerbgs[nwmers+i];
    for (j=0; j<w; j++) {
      probs[i]        += this->M[j][iseq[i+j]]       * this->mask[j];
      probs[nwmers+i] += this->M[w-j-1][3-iseq[i+j]] * this->mask[j];
      /* Complement of i = 3-i, convenient, eh? */
    }
    sumScores = log2_sum(sumScores,probs[i]);
    sumScores = log2_sum(sumScores,probs[nwmers+i]);
  }

  /* Renormalize according to Zoops model */
  float lambda      = gamma / (float(iseqlen)); /* Times 2 because of revcomplement */
  //printf("LAMBDA: %f\n",lambda);
  float log2lambda  = log(lambda) / log(2.0);
  float denominator = log2_sum( probebg + log(1.0-gamma)/log(2.0), sumScores + log(lambda)/log(2.0));
  //printf("C: probebg: %f + log(): %f, sumScores: %f + loglam() %f\n",
  //	 probebg, log(1.0-gamma)/log(2.0), sumScores, log2lambda);

  for (i=0; i<2*nwmers; i++) {
    //printf("C: Score: %f  + log_lambda: %f  - denom %f  = %f\n",
    //    probs[i],log2lambda,denominator,(probs[i] + log2lambda - denominator));
    Zs[i] = exp((probs[i] + log2lambda - denominator) *  0.693147180559945);
  }
  Zs[i] = -1e100; /* Terminator Kludge */
}

void SeqMat::zoops_m(int *iseq, int iseqlen, float *Zs, double newM[][_H]) {
  int i,j;
  int w           =  this->width;
  int nwmers      =  iseqlen - w + 1;

  for (i=0; i<nwmers; i++) {
    for (j=0; j<w; j++) {
      newM[  j  ][  iseq[i+j]] += Zs[i];
      newM[w-j-1][3-iseq[i+j]] += Zs[nwmers+i];
      /* Complement of i = 3-i, convenient, eh? */
    }
  }
}

void SeqMat::EMstep(Probelist *theProbelist, float gamma) {
  double newM[this->_W][this->_H];
  int i,j,jmax;
  float *probs =  new float[30000];
  float *Zs    =  new float[30000];
  float delta,diff;
  float sumZij;
  Probe *P;

  printf("#Mask: ");
  for (i=0;i<this->width;i++) printf("%3.1f ",this->mask[i]);
  printf("\n");
  
  int iterations = 0;
  for (iterations = 1;iterations<500; iterations++) {
    /* Zero New Matrix */
    for (i=0; i<this->width; i++) for (j=0; j<4; j++) newM[i][j]   = 0.0;
    for (i=0; i<this->width; i++) for (j=0; j<4; j++) this->counts[i][j] = 0.0;

    for (i=0; i < theProbelist->count; i++) {
      P = theProbelist->probes[i];
      //                                                              |-OUTPUT-|
      this->zoops_e(P->iseq,P->len,P->wmerbgs,P->probebg,this->gamma,  probs, Zs);
      this->zoops_m(P->iseq,P->len, Zs, this->counts);
      for (j=0; j<2*(P->len - this->width + 1); j++) {
	P->Zs[j] = Zs[j];
      }
    }

    /* Add pseudo counts and normalize */
    float coltot;
    float _t;
    for (i=0; i<this->width; i++) {
      coltot = 0;
      for (j=0; j<4; j++) coltot += this->counts[i][j];
      for (j=0; j<4; j++) {
	_t = this->counts[i][j];
	newM[i][j] = _t + (this->beta * this->bg[j] * coltot);         /* Add  */
	newM[i][j] = log(newM[i][j]/(coltot + this->beta)) / log(2.0); /* Norm */
      }
    }
    
    /* Compute the geometric distance from the last step */
    delta = 0;
    for (i=0; i<this->width; i++) for (j=0; j<4; j++) {
      diff  = pow(2.0,this->M[i][j]) - pow(2.0,newM[i][j]);
      delta += diff*diff;
      this->M[i][j] = newM[i][j];
    }
    delta = sqrt(delta);

    /* Re-estimate gamma */
    sumZij = 0;
    for (i=0; i < theProbelist->count; i++) {
      /* OLD CODE -- Ignored Zs of revcomp DBG 10-29-03    */
      /* for (j=0; j< theProbelist->probes[i]->len; j++) { */
      /*   sumZij += theProbelist->probes[i]->Zs[j];       */
      /* }                                                 */
      /* NEW CODE -- DBG 10-29-03                          */
      P = theProbelist->probes[i];
      jmax = 2*(P->len - this->width + 1);
      for (j=0; j<jmax; j++) {
	sumZij += P->Zs[j];
      }
    }

    /* this->gamma = variable gamma      */
    /* gamma       = fixed gamma (prior) */
    this->gamma = gamma_wt*gamma +  (1.0 - this->gamma_wt) * (sumZij/theProbelist->count);

    /* Print status */
    if (!(iterations%20)) printf("%3d: D: %f  G: %f\n",iterations,delta,this->gamma);
    if (!(iterations%10)) fflush(stdout);

    /* Have we converged? */
    if (delta <= this->deltamin) break;
  }

  this->joint = this->loglikelihood(theProbelist,this->gamma);

  delete[] probs;
  delete[] Zs;
}


float SeqMat::loglikelihood(Probelist *theProbelist, float gamma) {
  int i,j,k,w;
  float score, sum, Qi;
  float lambda;
  Probe *P;

  sum = 0;
  
  w = this->width;
  
  for (i=0; i<theProbelist->count; i++) {
    P = theProbelist->probes[i];
    lambda = gamma / P->len; // x2?? 
    Qi = 0;
    /* printf("%2d %3d\n",i,P->len-w+1); */
    for (j=0; j<P->len-w+1; j++) {
      score = 0;
      for (k=0; k<w; k++) {
	score += this->M[k][P->iseq[j+k]];        // * this->mask[k]? 
	score += this->M[w-k-1][3-P->iseq[j+k]];  // * this->mask[k]? 
      }
      sum = sum + P->Zs[j] + score;   // *2??  what about P->Zs[j+len?]
      Qi  = Qi + P->Zs[j];
    }
    sum = sum + (1-Qi) * P->probebg;
    sum = sum + (1-Qi) * log(1-gamma)/log(2.0);
    sum = sum + (  Qi) * log(lambda)/log(2.0);
    
  }
  return(sum);
}


  

