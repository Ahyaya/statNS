import numpy as np
import matplotlib.pyplot as plt

f4m = np.loadtxt("tovErr_rk4m_fixed.log", delimiter=',')
f5m = np.loadtxt("tovErr_rk5m_fixed.log", delimiter=',')
f5f = np.loadtxt("tovErr_rk5f_fixed.log", delimiter=',')
f5dp = np.loadtxt("tovErr_rk5dp_fixed.log", delimiter=',')
f5ck = np.loadtxt("tovErr_rk5ck_fixed.log", delimiter=',')

plt.figure()
ax=plt.gca()

ax.plot(f5ck[:,0], f5ck[:,1], "-", label = "Cash-Karp", color=plt.cm.Set2(0), linewidth=1, alpha=0.5)
ax.plot(f5f[:,0], f5f[:,1], "-", label = "RKF-5 (4)", color=plt.cm.Set2(1), linewidth=1, alpha=0.5)
ax.plot(f5m[:,0], f5m[:,1], "-", label = "Dormand-Prince (M)", color=plt.cm.Set2(2), linewidth=1, alpha=0.5)
ax.plot(f5dp[:,0], f5dp[:,1], "-", label = "Dormand-Prince", color=plt.cm.Set2(3), linewidth=1, alpha=0.5)
ax.plot(f4m[:,0], f4m[:,1], "k-", linewidth=0.5, label = "RK-4 (3)", alpha=0.5)

#ax.plot(np.array([0,11.5]), np.array([1e-8,1e-8]), "k--", linewidth = 2)
#ax.plot(np.array([0,11.5]), np.array([1e-10,1e-10]), "k--", linewidth = 2)

ax.minorticks_on()
ax.tick_params(which="both", top=1, right=1)
ax.tick_params(which="both", axis="both", direction="in", width=2.0)
ax.tick_params(which="major", length=6.0)
ax.tick_params(which="minor", length=4.0)

ax.spines["top"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)

ax.set_yscale("log")
ax.set_ylim(1e-13,1e-4)

labels = ax.get_xticklabels()+ax.get_yticklabels()
for label in labels:
    label.set_fontweight("bold")

ax.set_xlabel("Radius [$\mathrm{km}$]", fontdict={"weight":"bold","size":12})
ax.set_ylabel("Local Error", fontdict={"weight":"bold","size":12})

ax.legend(loc="upper center", frameon=0, fontsize=9.6, borderpad=1.2, ncol=1)
#plt.savefig("err_fix_h1.pdf")
plt.show()