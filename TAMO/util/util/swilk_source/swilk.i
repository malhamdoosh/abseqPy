// Copyright (2005) Whitehead Institute for Biomedical Research
// All Rights Reserved
// Author: David Benjamin Gordon

/* File : swilk.i */
%module  swilk
%include cpointer.i
%include carrays.i
%{
  /* Put headers and other declarations here */
%}

extern void swilkf(int &init, float *x, int &n, int &n1, int &n2, 
	           float *a, float &w, float &pw, int &ifault);

%pointer_class(int, intp);
%pointer_class(float, floatp);
%array_class(float, floatArray);

%pythoncode %{
def swilk(vect):
    init_p   = intp()   ; init_p.assign(0)
    n_p      = intp()   ; n_p.assign(len(vect))
    n1_p     = intp()   ; n1_p.assign(len(vect))
    n2_p     = intp()   ; n2_p.assign(len(vect)/2)
    w_p      = floatp() 
    pw_p     = floatp() 
    status_p = intp() 
    zero     = [0] * len(vect)

    vect.sort()
    x_array = floatArray(len(vect))
    zero_array = floatArray(len(vect))
    for i in range(len(vect)):
        x_array[i]=vect[i]
        zero_array[i]=0
        
    swilkf(init_p, x_array, n_p ,n1_p ,n2_p ,zero_array, w_p, pw_p, status_p)
    w  = w_p.value()
    pw = pw_p.value()
    
    status = status_p.value()

    statuses = ['ok ',                                       #  0
		'n1 < 3 ',                                   #  1
		'n > 5000 (non-fatal error) ',               #  2
		'n2 < n/2, so insufficient storage for a()', #  3
		'n1 > n or ((n1 < n) & (n < 20)) ',          #  4
		'proportion censored, (n-n1)/n, > 0.8 ',     #  5
		'effectively zero range (assuming data sorted) ', #  6
		'data not in increasing sort order ']        #  7

    if status != 0:
        print 'SHAPIRO-WILK TEST Error:', statuses[status]
        wp = 0
    return w,pw
%}


