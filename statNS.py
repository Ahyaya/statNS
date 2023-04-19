from ctypes import *

#dynamically load C library libstatNS.so as object libstatNS,
#please specify the path of the .so library according to your current directory
libstatNS=cdll.LoadLibrary('./lib/libstatNS.so')

#Begin of C2py
#=====================================================================
#set a map to translate variable types from C function to Python function
class CompactStar_t(Structure):
    _fields_ = (("P",c_double),("r",c_double),("Rho",c_double),("M",c_double),("Ma",c_double),("Mp",c_double),("I",c_double),("Ag00",c_double),("y",c_double),("k2",c_double),("Lambda",c_double),("Vs",c_double),("freq",c_double),("dampTime",c_double))
class EoS_t(Structure):
    _fields_ = (("length",c_int),("RhomaxSI",c_double),("RhominSI",c_double),("Rhomax",c_double),("Rhomin",c_double),("Pmax",c_double),("Pmin",c_double),("Mmax",c_double),("Rhoc_MmaxSI",c_double),("lgRho",c_double*2048),("lgP",c_double*2048),("lgRho_SI",c_double*2048),("lgP_SI",c_double*2048),("FilePath",c_char*128),("lgDen",c_double*2048),("lgDen_SI",c_double*2048),("XL",c_double),("Ksym",c_double),("Jsym",c_double),("J0",c_double))
class EoS_opt_t(Structure):
    _fields_ = (("secure",c_int),("causality",c_int),("total_length",c_int),("core_length",c_int),("dU",c_double),("maxU",c_double))

libstatNS.loadEoS.argtypes = (POINTER(EoS_t),c_char_p)
libstatNS.loadEoS.restype = c_int
loadEoS = libstatNS.loadEoS

libstatNS.fmode_mt.argtypes = (POINTER(CompactStar_t), POINTER(EoS_t), POINTER(c_double), c_int, c_int)
libstatNS.fmode_mt.restype = c_int
fmode_mt = libstatNS.fmode_mt

libstatNS.solveTOV_mt.argtypes = (POINTER(CompactStar_t), POINTER(EoS_t), POINTER(c_double), c_int, c_int)
libstatNS.solveTOV_mt.restype = c_int
solveTOV_mt = libstatNS.solveTOV_mt

libstatNS.solveTOV.argtypes = (POINTER(CompactStar_t), POINTER(EoS_t), c_double)
libstatNS.solveTOV.restype = c_int
solveTOV = libstatNS.solveTOV

libstatNS.getM_s.argtypes = (POINTER(EoS_t), c_double)
libstatNS.getM_s.restype = c_double
getM_s = libstatNS.getM_s
getM = getM_s

libstatNS.getM_fm.argtypes = (POINTER(EoS_t), c_double)
libstatNS.getM_fm.restype = c_double
getM_fm = libstatNS.getM_fm

libstatNS.getM_s_mt.argtypes = (POINTER(c_double), POINTER(EoS_t), POINTER(c_double), c_int, c_int)
libstatNS.getM_s_mt.restype = c_int
getM_s_mt = libstatNS.getM_s_mt
getM_mt = getM_s_mt

libstatNS.getM_fm_mt.argtypes = (POINTER(c_double), POINTER(EoS_t), POINTER(c_double), c_int, c_int)
libstatNS.getM_fm_mt.restype = c_int
getM_fm_mt = libstatNS.getM_fm_mt

libstatNS.getMmax.argtypes = (POINTER(EoS_t),)
libstatNS.getMmax.restype = c_double
getMmax = libstatNS.getMmax

libstatNS.getMmax_fm.argtypes = (POINTER(EoS_t),)
libstatNS.getMmax_fm.restype = c_double
getMmax_fm = libstatNS.getMmax_fm

libstatNS.M2Rhoc.argtypes = (POINTER(EoS_t), c_double)
libstatNS.M2Rhoc.restype = c_double
M2Rhoc = libstatNS.M2Rhoc

libstatNS.M2Rhoc_Arr_s.argtypes=(POINTER(c_double), POINTER(EoS_t), POINTER(c_double), c_int, c_int)
libstatNS.M2Rhoc_Arr_s.restype = c_int
M2Rhoc_Arr_s = libstatNS.M2Rhoc_Arr_s
M2Rhoc_Arr = M2Rhoc_Arr_s

libstatNS.M2Rhoc_Arr_fm.argtypes=(POINTER(c_double), POINTER(EoS_t), POINTER(c_double), c_int, c_int)
libstatNS.M2Rhoc_Arr_fm.restype = c_int
M2Rhoc_Arr_fm = libstatNS.M2Rhoc_Arr_fm

libstatNS.set_EoS_default_opt.argtypes = (POINTER(EoS_opt_t),)
libstatNS.set_EoS_default_opt.restype = c_int
set_EoS_default_opt = libstatNS.set_EoS_default_opt

libstatNS.genAmEoS.argtypes = (POINTER(EoS_t), c_double, c_double, c_double, c_double, POINTER(EoS_opt_t))
libstatNS.genAmEoS.restype = c_int
genAmEoS = libstatNS.genAmEoS

libstatNS.saveEoS.argtypes = (POINTER(EoS_t), c_char_p)
libstatNS.saveEoS.restype = c_int
saveEoS = libstatNS.saveEoS

#End of C2Py
#=====================================================================
