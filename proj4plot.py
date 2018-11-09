from pylab import *
rc('text', usetex=True)
rc('font', family='serif')

#
# filename = "4b2.txt"
# dat = loadtxt(filename, skiprows=1)
# anal = array([-1.99598, 0.998661, 0.0320823, 3.9933])
#
# E = dat[:, 1]
# M = dat[:, 5]
# Cv = dat[:, 7]
# chi = dat[:, 6]
# cyc = dat[:, -1]
#
#
# def absE(A, S):
#     # print A
#     # print S
#     return abs((A - S)/A)
#
#
# errE = zeros(len(E))
# errM = zeros(len(E))
# errCv = zeros(len(E))
# errChi = zeros(len(E))
#
# errE = absE(anal[0], E)
# errM = absE(anal[1], M)
# errCv = absE(anal[2], Cv)
# errChi = absE(anal[3], chi)
#
f_size = 28
#
# # legend(fontsize=f_size-2)
# #
# # xticks(size=f_size-2, rotation=30)
# # yticks(size=f_size-2, rotation=30)
# # ylabel(r'Absolute error', size=f_size)
# # xlabel(r'Number of Monte Carlo cycles', size=f_size)
# # tight_layout()
# #
# # show()

# loglog(cyc, errCv, 'r',label=r"Heat capacity $C_v$")
# loglog(cyc, errChi, 'b',label=r"Susceptibility $\chi$")
# loglog(cyc, errE, 'g',label=r"Energy $\langle E \rangle $")
# loglog(cyc, errM, 'k', label=r"Magnetization $|M|$")
# legend(fontsize=f_size-6)
# xticks(size=f_size-2, rotation=30)
# yticks(size=f_size-2, rotation=30)
# ylabel(r'Absolute error', size=f_size)
# xlabel(r'Number of Monte Carlo cycles', size=f_size)
# tight_layout()
# show()
###############################################################################

# cyc = 200
# cyc_string = str(cyc)
#
# ER_filename = "0calibrate20Cycles"+ cyc_string +"Ordered0.bin"
# EO_filename = "0calibrate20Cycles"+ cyc_string +"Ordered1.bin"
# MR_filename = "1calibrate20Cycles"+ cyc_string +"Ordered0.bin"
# MO_filename = "1calibrate20Cycles"+ cyc_string +"Ordered1.bin"
#
# cycles = np.linspace(1, cyc, cyc+2)
# ERvalues = np.fromfile(ER_filename, dtype=np.int32)/400.0
# EOvalues = np.fromfile(EO_filename, dtype=np.int32)/400.0
# MRvalues = np.fromfile(MR_filename, dtype=np.int32)/400.0
# MOvalues = np.fromfile(MO_filename, dtype=np.int32)/400.0

#
# plot(cycles, ERvalues,'r', label="Random matrix")
# plot(cycles, EOvalues,'b', label="Ordered matrix")
# # plt.plot(cycles, MRvalues)
# # plt.plot(cycles, MOvalues)
# xlabel(r"Number of Monte Carlo cycles", size=f_size)
# ylabel(r"$\langle E \rangle$", size=f_size)
# xticks(size=f_size-2, rotation=30)
# yticks(size=f_size-2, rotation=30)
# legend(fontsize=f_size-6)
# # ylim(-2.06,-1)
# # plt.axis([0, cyc+2, -2.1, -1.5 ])
# tight_layout()
# show()
#
# plot(cycles, MRvalues,'r', label="Random matrix")
# plot(cycles, MOvalues,'b', label="Ordered matrix")
# # plt.plot(cycles, MRvalues)
# # plt.plot(cycles, MOvalues)
# xlabel(r"Number of Monte Carlo cycles", size=f_size)
# ylabel(r"$\langle M \rangle$", size=f_size)
# xticks(size=f_size-2, rotation=30)
# yticks(size=f_size-2, rotation=30)
# legend(fontsize=f_size-6)
# # ylim(-2.06,-1)
# # plt.axis([0, cyc+2, -2.1, -1.5 ])
# tight_layout()
# show()
##############################################################################
f_size = 24
data = np.loadtxt("acceptsRatioVsT.txt")

T = data[:, 0]
ratio = data[:, 1]

plot(T, ratio)
xlabel(r"$T$ [kT/J]", size=f_size)
ylabel("Acceptance ratio", size=f_size)
xticks(size=f_size-2, rotation=30)
yticks(size=f_size-2, rotation=30)
tight_layout()
savefig("acceptRatio.pdf")
show()
