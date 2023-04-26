import numpy as np
import matplotlib.pyplot as plt
import statNS

myEoS = statNS.EoS_t()
Results = (statNS.CompactStar_t * 16)()
rhocArray = (statNS.c_double * 16)()

# Be careful the b"xxx" format is necessary in Python, strange isn't it?
statNS.loadEoS(statNS.pointer(myEoS), b"./EoS_lib/APR.txt")
Mmax = statNS.getMmax(statNS.pointer(myEoS))
Rhomin = statNS.M2Rhoc(statNS.pointer(myEoS), 1.05)
rhocArray = np.linspace(Rhomin, myEoS.Rhoc_MmaxSI, 16).ctypes.data_as(statNS.POINTER(statNS.c_double))
statNS.solveTOV_mt(Results, statNS.pointer(myEoS), rhocArray, 16, 4)

APR_R = np.zeros(16)
APR_M = np.zeros(16)
for pf in range(0, 16):
    APR_R[pf] = Results[pf].r
    APR_M[pf] = Results[pf].M

statNS.loadEoS(statNS.pointer(myEoS), b"./EoS_lib/ALF2.txt")
Mmax = statNS.getMmax(statNS.pointer(myEoS))
Rhomin = statNS.M2Rhoc(statNS.pointer(myEoS), 1.05)
rhocArray = np.linspace(Rhomin, myEoS.Rhoc_MmaxSI, 16).ctypes.data_as(statNS.POINTER(statNS.c_double))
statNS.solveTOV_mt(Results, statNS.pointer(myEoS), rhocArray, 16, 4)

ALF_R = np.zeros(16)
ALF_M = np.zeros(16)
for pf in range(0, 16):
    ALF_R[pf] = Results[pf].r
    ALF_M[pf] = Results[pf].M

statNS.loadEoS(statNS.pointer(myEoS), b"./EoS_lib/MPA1.txt")
Mmax = statNS.getMmax(statNS.pointer(myEoS))
Rhomin = statNS.M2Rhoc(statNS.pointer(myEoS), 1.05)
rhocArray = np.linspace(Rhomin, myEoS.Rhoc_MmaxSI, 16).ctypes.data_as(statNS.POINTER(statNS.c_double))
statNS.solveTOV_mt(Results, statNS.pointer(myEoS), rhocArray, 16, 4)

MPA_R = np.zeros(16)
MPA_M = np.zeros(16)
for pf in range(0, 16):
    MPA_R[pf] = Results[pf].r
    MPA_M[pf] = Results[pf].M

plt.figure()
ax=plt.gca()

ax.plot(APR_R, APR_M, "-s", label = "APR", color=plt.cm.Set2(0), linewidth=2, alpha=0.5)
ax.plot(ALF_R, ALF_M, "-^", label = "MPA1", color=plt.cm.Set2(1), linewidth=2, alpha=0.5)
ax.plot(MPA_R, MPA_M, "-o", label = "ALF2", color=plt.cm.Set2(2), linewidth=2, alpha=0.5)

ax.minorticks_on()
ax.tick_params(which="both", top=1, right=1)
ax.tick_params(which="both", axis="both", direction="in", width=2.0)
ax.tick_params(which="major", length=6.0)
ax.tick_params(which="minor", length=4.0)

ax.spines["top"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)

labels = ax.get_xticklabels()+ax.get_yticklabels()
for label in labels:
    label.set_fontweight("bold")

ax.set_xlim(9.45,14.5)
ax.set_ylim(0.86,2.65)

ax.set_xlabel("Radius [$\mathrm{km}$]", fontdict={"weight":"bold","size":12})
ax.set_ylabel("Mass [$\mathrm{M_\odot}$]", fontdict={"weight":"bold","size":12})

ax.legend(loc="upper right", frameon=0, fontsize=9.6, borderpad=1.2, ncol=1)
#plt.savefig("RM_example.pdf")
plt.show()