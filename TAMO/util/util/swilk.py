# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.
import _swilk
def _swig_setattr(self,class_type,name,value):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    self.__dict__[name] = value

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


swilkf = _swilk.swilkf

class intp(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, intp, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, intp, name)
    def __init__(self,*args):
        _swig_setattr(self, intp, 'this', apply(_swilk.new_intp,args))
        _swig_setattr(self, intp, 'thisown', 1)
    def __del__(self, destroy= _swilk.delete_intp):
        try:
            if self.thisown: destroy(self)
        except: pass
    def assign(*args): return apply(_swilk.intp_assign,args)
    def value(*args): return apply(_swilk.intp_value,args)
    def cast(*args): return apply(_swilk.intp_cast,args)
    __swig_getmethods__["frompointer"] = lambda x: _swilk.intp_frompointer
    if _newclass:frompointer = staticmethod(_swilk.intp_frompointer)
    def __repr__(self):
        return "<C intp instance at %s>" % (self.this,)

class intpPtr(intp):
    def __init__(self,this):
        _swig_setattr(self, intp, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, intp, 'thisown', 0)
        _swig_setattr(self, intp,self.__class__,intp)
_swilk.intp_swigregister(intpPtr)
intp_frompointer = _swilk.intp_frompointer


class floatp(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, floatp, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, floatp, name)
    def __init__(self,*args):
        _swig_setattr(self, floatp, 'this', apply(_swilk.new_floatp,args))
        _swig_setattr(self, floatp, 'thisown', 1)
    def __del__(self, destroy= _swilk.delete_floatp):
        try:
            if self.thisown: destroy(self)
        except: pass
    def assign(*args): return apply(_swilk.floatp_assign,args)
    def value(*args): return apply(_swilk.floatp_value,args)
    def cast(*args): return apply(_swilk.floatp_cast,args)
    __swig_getmethods__["frompointer"] = lambda x: _swilk.floatp_frompointer
    if _newclass:frompointer = staticmethod(_swilk.floatp_frompointer)
    def __repr__(self):
        return "<C floatp instance at %s>" % (self.this,)

class floatpPtr(floatp):
    def __init__(self,this):
        _swig_setattr(self, floatp, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, floatp, 'thisown', 0)
        _swig_setattr(self, floatp,self.__class__,floatp)
_swilk.floatp_swigregister(floatpPtr)
floatp_frompointer = _swilk.floatp_frompointer


class floatArray(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, floatArray, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, floatArray, name)
    def __init__(self,*args):
        _swig_setattr(self, floatArray, 'this', apply(_swilk.new_floatArray,args))
        _swig_setattr(self, floatArray, 'thisown', 1)
    def __del__(self, destroy= _swilk.delete_floatArray):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __getitem__(*args): return apply(_swilk.floatArray___getitem__,args)
    def __setitem__(*args): return apply(_swilk.floatArray___setitem__,args)
    def cast(*args): return apply(_swilk.floatArray_cast,args)
    __swig_getmethods__["frompointer"] = lambda x: _swilk.floatArray_frompointer
    if _newclass:frompointer = staticmethod(_swilk.floatArray_frompointer)
    def __repr__(self):
        return "<C floatArray instance at %s>" % (self.this,)

class floatArrayPtr(floatArray):
    def __init__(self,this):
        _swig_setattr(self, floatArray, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, floatArray, 'thisown', 0)
        _swig_setattr(self, floatArray,self.__class__,floatArray)
_swilk.floatArray_swigregister(floatArrayPtr)
floatArray_frompointer = _swilk.floatArray_frompointer


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



