// Copyright (2005) Whitehead Institute for Biomedical Research
// All Rights Reserved
// Author: David Benjamin Gordon

#ifdef __cplusplus
extern "C" {
#endif
#include "Python.h"
#ifdef __cplusplus
}
#endif

#include <vector>
#include <sys/types.h>
#include <regex.h>
#include <stdlib.h>
#include <stdio.h>

typedef float*  flist;
typedef double* dlist;

int    *seq2int    (const char *seq, int len = 0);
void    print_seq  (int *iseq, int len);
double *list2double(dlist dptr);
void    printdouble(double *dlist, int len);
float   log2_sum   (float logx,float logy);

class Probe {
 public:
  int *iseq;
  int len;
  float probebg;
  double *wmerbgs;
  float *Zs;

 public:
  Probe(int *_iseq, int _len, float _probebg, double *_wmerbgs) {
    this->iseq    = _iseq;
    this->len     = _len;
    this->probebg = _probebg;
    this->wmerbgs = _wmerbgs;
    this->Zs      =  new float[_len*2+1];
  }
  ~Probe() {
    delete(this->Zs);
  }
};


class Probelist_str {
 private:
 public:
  int count;
  char *probes[22000];

 public:
  Probelist_str() {this->count = 0;} 
  void append(char *seq) {
    int len;
    for (len=0; seq[len] != 0L; len++);
    probes[this->count] = new char[len+10];
    strcpy(probes[this->count],seq);
    this->count++;
  }
  ~Probelist_str() {
    for (int i=0; i<this->count; i++) {
      delete this->probes[i];
    }
  }
  int count_re_matches(char *expr, int len=0) {
    int i;
    int n = 0;
    int result;
    regex_t *re_compiled = NULL;
    //regex_t *re_compiled = (regex_t *) malloc(sizeof(regex_t));
    //regmatch_t pmatches[10];
    regcomp(re_compiled,expr,REG_NOSUB);
    for (i=0; i<this->count; i++) {
      result = regexec(re_compiled,this->probes[i],0,NULL,0);
      if (result == 0) n++;
    }
    regfree(re_compiled);
    //free(re_compiled);
    return n;
  }
};

class Probelist {
 private:
  
 public:
  int count;
  Probe *probes[22000];

 public:
  Probelist() {this->count = 0;};
  void append(int *iseq,int len, float probebg, double *wmerbgs) {
    probes[this->count] = new Probe(iseq,len,probebg,wmerbgs);
    this->count++;
  }
  ~Probelist() {
    for (int i=0; i<this->count; i++) {
      delete this->probes[i];
    }
  }

  float get_Z(int i, int j) { 
    float ans;
    if (i < this->count && j < 2*(this->probes[i]->len)) { /*Approximate bounds */
      ans = this->probes[i]->Zs[j];
    } else {
      printf ("ERROR: Request for Zij from probe %d, element %d, (>%d or >%d)\n",
	      i,j,this->count,this->probes[i]->len);
      ans = 0;
    }
    return(ans);
  }

  flist get_Zlist(int probenum) {  /* I don't handle memory correctly.  Don't use! */
    flist ans;
    if (probenum < this->count) {
      ans = this->probes[probenum]->Zs;
    } else {
      printf ("Request for Zlist from probe %d, which does not exist (>%d)\n",
	      probenum,this->count);
      ans = NULL;
    }
    return(ans);
  }

};


class SeqMat {
 private:
  enum   {_W = 200, _H = 5};
  double M[_W][_H];
  double counts[_W][_H];

 public:  //Variables
  int    width;
  float  gamma;
  float  gamma_wt;
  float  deltamin;
  float  beta;
  float  bg[4];
  float  mask[100];
  double joint;

 public: //Functions
  SeqMat(int w) {_init_(w);};
  double scanbest(const char *seq);
  double sumscoresabove(char *seq, float cutoff);
  double score(int *idxarray);
  void   set  (int i, int j, double value);
  double get  (int i, int j);
  double get_c(int i, int j);
  void   compute_ambig();
  std::vector<int> matchstarts(char *seq, float cutoff);
  flist  score_probe(char *seq);
  void EMstep(Probelist *theProbelist, float gamma);
  void setBg(double A, double C, double G, double T);
  void setmask(int pos, float f);
  float loglikelihood(Probelist *theProbelist, float gamma);

  /*
  PyObject* ZoopsE_probe(int *iseq, int iseqlen,
			 double *wmerbgs, 
			 float gamma, float probebg);
  */



 private:
  void zoops_e(int *iseq, int iseqlen, double *wmerbgs, 
	       float probebg, float gamma,
	       float *probs, float *Zs);
  void zoops_m(int *iseq, int iseqlen, float *Zs, double newM[][_H]);


 private: //Functions
  void  _init_(int w);
};

