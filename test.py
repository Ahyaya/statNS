from ctypes import *
import numpy as np

#dynamically load C library libstatNS.so as object statNS,
statNS=cdll.LoadLibrary('./lib/libstatNS.so')

#Begin of C2py
#set a map to translate variable types from C function to Python function
class CompactStar_t(Structure):
    _fields_ = (("P",c_double),("r",c_double),("Rho",c_double),("M",c_double),("Ma",c_double),("Mp",c_double),("I",c_double),("Ag00",c_double),("y",c_double),("k2",c_double),("Lambda",c_double),("Vs",c_double),("freq",c_double),("dampTime",c_double))
class EoS_t(Structure):
    _fields_ = (("length",c_int),("RhomaxSI",c_double),("RhominSI",c_double),("Rhomax",c_double),("Rhomin",c_double),("Pmax",c_double),("Pmin",c_double),("Mmax",c_double),("Rhoc_MmaxSI",c_double),("lgRho",c_double*3000),("lgP",c_double*3000),("lgRho_SI",c_double*3000),("lgP_SI",c_double*3000),("FilePath",c_char*128))
statNS.loadEoS.argtypes = (POINTER(EoS_t),c_char_p)
statNS.loadEoS.restype = c_int
statNS.fmode_mt.argtypes = (POINTER(CompactStar_t), POINTER(EoS_t), POINTER(c_double), c_int, c_int)
statNS.fmode_mt.restype = c_int
statNS.solveTOV_mt.argtypes = (POINTER(CompactStar_t), POINTER(EoS_t), POINTER(c_double), c_int, c_int)
statNS.solveTOV_mt.restype = c_int
statNS.solveTOV.argtypes = (POINTER(CompactStar_t), POINTER(EoS_t), c_double)
statNS.solveTOV.restype = c_int
statNS.M2Rhoc_Arr_fm.argtypes=(POINTER(c_double), POINTER(EoS_t), POINTER(c_double), c_int, c_int);
statNS.M2Rhoc_Arr_fm.restype = c_int
statNS.M2Rhoc_Arr_s.argtypes=(POINTER(c_double), POINTER(EoS_t), POINTER(c_double), c_int, c_int);
statNS.M2Rhoc_Arr_s.restype = c_int
statNS.M2Rhoc.argtypes = (POINTER(EoS_t), c_double)
statNS.M2Rhoc.restype = c_double
#End of C2Py

#=====================================================================
#Example 0: Compute f-mode frequencies given an array of central density

#define central density sequence in the SI unit kg/m^3
RhocSI = (c_double*19)(6.5e17,7e17,7.2e17,8e17,8.2e17,8.4e17,8.6e17,8.8e17,9e17,9.2e17,9.5e17,9.8e17,1e18,1.05e18,1.08e18,1.1e18,1.12e18,1.14e18,1.2e18)

#initiate an EoS_t type variable to store your EoS info
myEoS = EoS_t()

#initiate a CompactStar_t type array for storing the output results,
#you may register a large enough array to store results when using longer input.
Results = (CompactStar_t * 18)()

#load EOS file "APR.txt" to myEoS, be careful the b"xxx" format is Python dedicated, strange isn't it?
statNS.loadEoS(pointer(myEoS), b"APR.txt")

#call f-mode multi-thread computation, it will return the results to Results[],
#the last two numbers represent: 18 the length of RhocSI[], 5 the total threads for computation.
statNS.fmode_mt(Results,pointer(myEoS),RhocSI,18,5)

#print results to console
print("\nExample 0 output:")
for pf in range(0,18):
    print("Rhoc={:.5e}, M={:.2f}, R={:.2f}, freq={:.2f}, dmpTime={:.4f}".format(Results[pf].Rho,Results[pf].M,Results[pf].r,Results[pf].freq,Results[pf].dampTime))

#=====================================================================
#Example 1: Compute f-mode frequencies given an array of mass

#define mass array sequence in the unit of solar mass
massArray = (c_double*16)(1.25, 1.325, 1.4, 1.44 ,1.46, 1.48, 1.5, 1.55, 1.62, 1.64, 1.67, 1.7123, 1.7456, 1.8257, 1.90123, 1.95)

#This rhocArray[] is used to store the central density, it needs a proper size not less than massArray[]
rhocArray = (c_double*16)()

#This function compute the central density that corresponds to the mass array one by one,
#The central density will be written into rhocArray[],
#the last two numbers are, 16 the length fo massArray[], 8 the desired threads to use
statNS.M2Rhoc_Arr_fm(rhocArray, pointer(myEoS), massArray, 16, 8);

#Call the f-mode multi-thread computation with rhocArray[], see also in Example 0.
statNS.fmode_mt(Results, pointer(myEoS), rhocArray, 16, 8);

#print results to console
print("\nExample 1 output:")
for pf in range(0,16):
    print("Rhoc={:.5e} -->\tM={:.4f} freq={:.4f} dTime={:.4f}".format(rhocArray[pf],Results[pf].M,Results[pf].freq,Results[pf].dampTime))


print("Done\n")
