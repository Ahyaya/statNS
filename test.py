from ctypes import *
import numpy as np

statNS=cdll.LoadLibrary('./lib/libstatNS.so')

class CompactStar_t(Structure):
    _fields_ = (("P",c_double),("r",c_double),("Rho",c_double),("M",c_double),("Ma",c_double),("Mp",c_double),("I",c_double),("Ag00",c_double),("y",c_double),("k2",c_double),("Lambda",c_double),("Vs",c_double),("freq",c_double),("dampTime",c_double))

statNS.fmode_mt.argtypes = (c_int, POINTER(CompactStar_t), POINTER(c_double), c_int)
statNS.fmode_mt.restype = c_int

statNS.solveTOV_mt.argtypes = (c_int, POINTER(CompactStar_t), POINTER(c_double), c_int)
statNS.solveTOV_mt.restype = c_int

statNS.solveTOV.argtypes = (c_double, POINTER(CompactStar_t))
statNS.solveTOV.restype = c_double

statNS.loadEoS.argtypes = (c_char_p,)
statNS.loadEoS.restype = c_int

statNS.M2Rhoc.argtypes = (c_double,)
statNS.M2Rhoc.restype = c_double

cdenArray = (c_double*19)(6.5e17,7e17,7.2e17,8e17,8.2e17,8.4e17,8.6e17,8.8e17,9e17,9.2e17,9.5e17,9.8e17,1e18,1.05e18,1.08e18,1.1e18,1.12e18,1.14e18,1.2e18)

Results = (CompactStar_t * 18)()

statNS.loadEoS(b"APR.txt")
statNS.fmode_mt(4,Results,cdenArray,18)
print("fmode:\n")
for pf in range(0,18):
    print("rho=",Results[pf].Rho,"M=",Results[pf].M,"freq=",Results[pf].freq)

statNS.solveTOV_mt(4,Results,cdenArray,18)
print("\nsolveTOV:\n")
for pf in range(0,18):
    print("rho=",Results[pf].Rho,"M=",Results[pf].M,"I=",Results[pf].I,"Lambda=",Results[pf].Lambda)

print("final test: searching for M=1.44007")
myRhoc = statNS.M2Rhoc(1.44007)
statNS.solveTOV(myRhoc,Results)
print("M=", Results[0].M, "R=", Results[0].r, "I=", Results[0].I, "Lambda=", Results[0].Lambda)
