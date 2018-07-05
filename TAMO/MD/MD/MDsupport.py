# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.
import _MDsupport
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


class IntVector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, IntVector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, IntVector, name)
    def __init__(self,*args):
        _swig_setattr(self, IntVector, 'this', apply(_MDsupport.new_IntVector,args))
        _swig_setattr(self, IntVector, 'thisown', 1)
    def __len__(*args): return apply(_MDsupport.IntVector___len__,args)
    def __nonzero__(*args): return apply(_MDsupport.IntVector___nonzero__,args)
    def clear(*args): return apply(_MDsupport.IntVector_clear,args)
    def append(*args): return apply(_MDsupport.IntVector_append,args)
    def pop(*args): return apply(_MDsupport.IntVector_pop,args)
    def __getitem__(*args): return apply(_MDsupport.IntVector___getitem__,args)
    def __getslice__(*args): return apply(_MDsupport.IntVector___getslice__,args)
    def __setitem__(*args): return apply(_MDsupport.IntVector___setitem__,args)
    def __setslice__(*args): return apply(_MDsupport.IntVector___setslice__,args)
    def __delitem__(*args): return apply(_MDsupport.IntVector___delitem__,args)
    def __delslice__(*args): return apply(_MDsupport.IntVector___delslice__,args)
    def __del__(self, destroy= _MDsupport.delete_IntVector):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __repr__(self):
        return "<C IntVector instance at %s>" % (self.this,)

class IntVectorPtr(IntVector):
    def __init__(self,this):
        _swig_setattr(self, IntVector, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, IntVector, 'thisown', 0)
        _swig_setattr(self, IntVector,self.__class__,IntVector)
_MDsupport.IntVector_swigregister(IntVectorPtr)

class DoubleVector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DoubleVector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DoubleVector, name)
    def __init__(self,*args):
        _swig_setattr(self, DoubleVector, 'this', apply(_MDsupport.new_DoubleVector,args))
        _swig_setattr(self, DoubleVector, 'thisown', 1)
    def __len__(*args): return apply(_MDsupport.DoubleVector___len__,args)
    def __nonzero__(*args): return apply(_MDsupport.DoubleVector___nonzero__,args)
    def clear(*args): return apply(_MDsupport.DoubleVector_clear,args)
    def append(*args): return apply(_MDsupport.DoubleVector_append,args)
    def pop(*args): return apply(_MDsupport.DoubleVector_pop,args)
    def __getitem__(*args): return apply(_MDsupport.DoubleVector___getitem__,args)
    def __getslice__(*args): return apply(_MDsupport.DoubleVector___getslice__,args)
    def __setitem__(*args): return apply(_MDsupport.DoubleVector___setitem__,args)
    def __setslice__(*args): return apply(_MDsupport.DoubleVector___setslice__,args)
    def __delitem__(*args): return apply(_MDsupport.DoubleVector___delitem__,args)
    def __delslice__(*args): return apply(_MDsupport.DoubleVector___delslice__,args)
    def __del__(self, destroy= _MDsupport.delete_DoubleVector):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __repr__(self):
        return "<C DoubleVector instance at %s>" % (self.this,)

class DoubleVectorPtr(DoubleVector):
    def __init__(self,this):
        _swig_setattr(self, DoubleVector, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, DoubleVector, 'thisown', 0)
        _swig_setattr(self, DoubleVector,self.__class__,DoubleVector)
_MDsupport.DoubleVector_swigregister(DoubleVectorPtr)

def Motif2c_PSSM(motif):
    c_PSSM = SeqMat(motif.width)
    LjT = zip(['A','C','G','T'],[0,1,2,3])
    for i in range(motif.width):
        for L,j in LjT:
            c_PSSM.set(i,j, motif.ll[i][L])
    c_PSSM.compute_ambig()
    return c_PSSM


seq2int = _MDsupport.seq2int

print_seq = _MDsupport.print_seq

list2double = _MDsupport.list2double

printdouble = _MDsupport.printdouble

log2_sum = _MDsupport.log2_sum

class Probe(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Probe, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Probe, name)
    __swig_setmethods__["iseq"] = _MDsupport.Probe_iseq_set
    __swig_getmethods__["iseq"] = _MDsupport.Probe_iseq_get
    if _newclass:iseq = property(_MDsupport.Probe_iseq_get,_MDsupport.Probe_iseq_set)
    __swig_setmethods__["len"] = _MDsupport.Probe_len_set
    __swig_getmethods__["len"] = _MDsupport.Probe_len_get
    if _newclass:len = property(_MDsupport.Probe_len_get,_MDsupport.Probe_len_set)
    __swig_setmethods__["probebg"] = _MDsupport.Probe_probebg_set
    __swig_getmethods__["probebg"] = _MDsupport.Probe_probebg_get
    if _newclass:probebg = property(_MDsupport.Probe_probebg_get,_MDsupport.Probe_probebg_set)
    __swig_setmethods__["wmerbgs"] = _MDsupport.Probe_wmerbgs_set
    __swig_getmethods__["wmerbgs"] = _MDsupport.Probe_wmerbgs_get
    if _newclass:wmerbgs = property(_MDsupport.Probe_wmerbgs_get,_MDsupport.Probe_wmerbgs_set)
    __swig_setmethods__["Zs"] = _MDsupport.Probe_Zs_set
    __swig_getmethods__["Zs"] = _MDsupport.Probe_Zs_get
    if _newclass:Zs = property(_MDsupport.Probe_Zs_get,_MDsupport.Probe_Zs_set)
    def __init__(self,*args):
        _swig_setattr(self, Probe, 'this', apply(_MDsupport.new_Probe,args))
        _swig_setattr(self, Probe, 'thisown', 1)
    def __del__(self, destroy= _MDsupport.delete_Probe):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __repr__(self):
        return "<C Probe instance at %s>" % (self.this,)

class ProbePtr(Probe):
    def __init__(self,this):
        _swig_setattr(self, Probe, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Probe, 'thisown', 0)
        _swig_setattr(self, Probe,self.__class__,Probe)
_MDsupport.Probe_swigregister(ProbePtr)

class Probelist_str(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Probelist_str, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Probelist_str, name)
    __swig_setmethods__["count"] = _MDsupport.Probelist_str_count_set
    __swig_getmethods__["count"] = _MDsupport.Probelist_str_count_get
    if _newclass:count = property(_MDsupport.Probelist_str_count_get,_MDsupport.Probelist_str_count_set)
    __swig_setmethods__["probes"] = _MDsupport.Probelist_str_probes_set
    __swig_getmethods__["probes"] = _MDsupport.Probelist_str_probes_get
    if _newclass:probes = property(_MDsupport.Probelist_str_probes_get,_MDsupport.Probelist_str_probes_set)
    def __init__(self,*args):
        _swig_setattr(self, Probelist_str, 'this', apply(_MDsupport.new_Probelist_str,args))
        _swig_setattr(self, Probelist_str, 'thisown', 1)
    def append(*args): return apply(_MDsupport.Probelist_str_append,args)
    def __del__(self, destroy= _MDsupport.delete_Probelist_str):
        try:
            if self.thisown: destroy(self)
        except: pass
    def count_re_matches(*args): return apply(_MDsupport.Probelist_str_count_re_matches,args)
    def __repr__(self):
        return "<C Probelist_str instance at %s>" % (self.this,)

class Probelist_strPtr(Probelist_str):
    def __init__(self,this):
        _swig_setattr(self, Probelist_str, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Probelist_str, 'thisown', 0)
        _swig_setattr(self, Probelist_str,self.__class__,Probelist_str)
_MDsupport.Probelist_str_swigregister(Probelist_strPtr)

class Probelist(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Probelist, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Probelist, name)
    __swig_setmethods__["count"] = _MDsupport.Probelist_count_set
    __swig_getmethods__["count"] = _MDsupport.Probelist_count_get
    if _newclass:count = property(_MDsupport.Probelist_count_get,_MDsupport.Probelist_count_set)
    __swig_setmethods__["probes"] = _MDsupport.Probelist_probes_set
    __swig_getmethods__["probes"] = _MDsupport.Probelist_probes_get
    if _newclass:probes = property(_MDsupport.Probelist_probes_get,_MDsupport.Probelist_probes_set)
    def __init__(self,*args):
        _swig_setattr(self, Probelist, 'this', apply(_MDsupport.new_Probelist,args))
        _swig_setattr(self, Probelist, 'thisown', 1)
    def append(*args): return apply(_MDsupport.Probelist_append,args)
    def __del__(self, destroy= _MDsupport.delete_Probelist):
        try:
            if self.thisown: destroy(self)
        except: pass
    def get_Z(*args): return apply(_MDsupport.Probelist_get_Z,args)
    def get_Zlist(*args): return apply(_MDsupport.Probelist_get_Zlist,args)
    def __repr__(self):
        return "<C Probelist instance at %s>" % (self.this,)

class ProbelistPtr(Probelist):
    def __init__(self,this):
        _swig_setattr(self, Probelist, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Probelist, 'thisown', 0)
        _swig_setattr(self, Probelist,self.__class__,Probelist)
_MDsupport.Probelist_swigregister(ProbelistPtr)

class SeqMat(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SeqMat, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SeqMat, name)
    __swig_setmethods__["width"] = _MDsupport.SeqMat_width_set
    __swig_getmethods__["width"] = _MDsupport.SeqMat_width_get
    if _newclass:width = property(_MDsupport.SeqMat_width_get,_MDsupport.SeqMat_width_set)
    __swig_setmethods__["gamma"] = _MDsupport.SeqMat_gamma_set
    __swig_getmethods__["gamma"] = _MDsupport.SeqMat_gamma_get
    if _newclass:gamma = property(_MDsupport.SeqMat_gamma_get,_MDsupport.SeqMat_gamma_set)
    __swig_setmethods__["gamma_wt"] = _MDsupport.SeqMat_gamma_wt_set
    __swig_getmethods__["gamma_wt"] = _MDsupport.SeqMat_gamma_wt_get
    if _newclass:gamma_wt = property(_MDsupport.SeqMat_gamma_wt_get,_MDsupport.SeqMat_gamma_wt_set)
    __swig_setmethods__["deltamin"] = _MDsupport.SeqMat_deltamin_set
    __swig_getmethods__["deltamin"] = _MDsupport.SeqMat_deltamin_get
    if _newclass:deltamin = property(_MDsupport.SeqMat_deltamin_get,_MDsupport.SeqMat_deltamin_set)
    __swig_setmethods__["beta"] = _MDsupport.SeqMat_beta_set
    __swig_getmethods__["beta"] = _MDsupport.SeqMat_beta_get
    if _newclass:beta = property(_MDsupport.SeqMat_beta_get,_MDsupport.SeqMat_beta_set)
    __swig_setmethods__["bg"] = _MDsupport.SeqMat_bg_set
    __swig_getmethods__["bg"] = _MDsupport.SeqMat_bg_get
    if _newclass:bg = property(_MDsupport.SeqMat_bg_get,_MDsupport.SeqMat_bg_set)
    __swig_setmethods__["mask"] = _MDsupport.SeqMat_mask_set
    __swig_getmethods__["mask"] = _MDsupport.SeqMat_mask_get
    if _newclass:mask = property(_MDsupport.SeqMat_mask_get,_MDsupport.SeqMat_mask_set)
    __swig_setmethods__["joint"] = _MDsupport.SeqMat_joint_set
    __swig_getmethods__["joint"] = _MDsupport.SeqMat_joint_get
    if _newclass:joint = property(_MDsupport.SeqMat_joint_get,_MDsupport.SeqMat_joint_set)
    def __init__(self,*args):
        _swig_setattr(self, SeqMat, 'this', apply(_MDsupport.new_SeqMat,args))
        _swig_setattr(self, SeqMat, 'thisown', 1)
    def scanbest(*args): return apply(_MDsupport.SeqMat_scanbest,args)
    def sumscoresabove(*args): return apply(_MDsupport.SeqMat_sumscoresabove,args)
    def score(*args): return apply(_MDsupport.SeqMat_score,args)
    def set(*args): return apply(_MDsupport.SeqMat_set,args)
    def get(*args): return apply(_MDsupport.SeqMat_get,args)
    def get_c(*args): return apply(_MDsupport.SeqMat_get_c,args)
    def compute_ambig(*args): return apply(_MDsupport.SeqMat_compute_ambig,args)
    def matchstarts(*args): return apply(_MDsupport.SeqMat_matchstarts,args)
    def score_probe(*args): return apply(_MDsupport.SeqMat_score_probe,args)
    def EMstep(*args): return apply(_MDsupport.SeqMat_EMstep,args)
    def setBg(*args): return apply(_MDsupport.SeqMat_setBg,args)
    def setmask(*args): return apply(_MDsupport.SeqMat_setmask,args)
    def loglikelihood(*args): return apply(_MDsupport.SeqMat_loglikelihood,args)
    def __del__(self, destroy= _MDsupport.delete_SeqMat):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __repr__(self):
        return "<C SeqMat instance at %s>" % (self.this,)

class SeqMatPtr(SeqMat):
    def __init__(self,this):
        _swig_setattr(self, SeqMat, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, SeqMat, 'thisown', 0)
        _swig_setattr(self, SeqMat,self.__class__,SeqMat)
_MDsupport.SeqMat_swigregister(SeqMatPtr)


