from pylab import *
rc('text', usetex=True)
rc('font', family='serif')

import matplotlib.mlab as mlab
from scipy.stats import norm

# # f_size = 22 [scale=0.55] for centered plots, 28 for neighboring plots [scale=0.44]
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
# f_size = 22
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
#
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
# ###########################################################################
# f_size = 28
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
#
# plot(cycles, ERvalues,'r', label="Random matrix")
# plot(cycles, EOvalues,'b', label="Ordered matrix")
# # plt.plot(cycles, MRvalues)
# # plt.plot(cycles, MOvalues)
# xlabel(r"Number of Monte Carlo cycles", size=f_size)
# ylabel(r"Energy", size=f_size)
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
# ylabel(r"Magnetization", size=f_size)
# xticks(size=f_size-2, rotation=30)
# yticks(size=f_size-2, rotation=30)
# legend(fontsize=f_size-6)
# # ylim(-2.06,-1)
# # plt.axis([0, cyc+2, -2.1, -1.5 ])
# tight_layout()
# show()
##########################################################################
# f_size = 26
# data = np.loadtxt("acceptsRatioVsT.txt")
#
# T = data[:, 0]
# ratio = data[:, 1]
#
# plot(T, ratio)
# xlabel(r" Temperature [$kT/J$]", size=f_size)
# ylabel("Number of accepted \n configurations", size=f_size)
# xticks(size=f_size-2, rotation=30)
# yticks(size=f_size-2, rotation=30)
# tight_layout()
# savefig("acceptRatio.pdf")
# show()
# ############################################################################


#############################################################################
f_size = 28
l_size = f_size - 9
# data = np.loadtxt("4e_dim40_cycles100000.txt", skiprows=1)
# data1 = np.loadtxt("4e_dim60_cycles100000.txt", skiprows=1)
# data2 = np.loadtxt("4e_dim80_cycles100000.txt", skiprows=1)
# data3 = np.loadtxt("4e_dim100_cycles100000.txt", skiprows=1)

# data = np.loadtxt("A_TEST_4e_dim40_cycles1000.txt", skiprows=1)
# data1 = np.loadtxt("A_TEST_4e_dim60_cycles1000.txt", skiprows=1)
# data2 = np.loadtxt("A_TEST_4e_dim80_cycles1000.txt", skiprows=1)
# data3 = np.loadtxt("A_TEST_4e_dim100_cycles1000.txt", skiprows=1)
#
# data = np.loadtxt("MONEY_4e_dim40_cycles100000.txt", skiprows=1)
# data1 = np.loadtxt("MONEY_4e_dim60_cycles100000.txt", skiprows=1)
# data2 = np.loadtxt("MONEY_4e_dim80_cycles100000.txt", skiprows=1)
# data3 = np.loadtxt("MONEY_4e_dim100_cycles100000.txt", skiprows=1)


data = np.loadtxt("NarrowT_4e_dim40_cycles100000.txt", skiprows=1)
data1 = np.loadtxt("NarrowT_4e_dim60_cycles100000.txt", skiprows=1)
data2 = np.loadtxt("NarrowT_4e_dim80_cycles100000.txt", skiprows=1)
data3 = np.loadtxt("NarrowT_4e_dim100_cycles100000.txt", skiprows=1)

T = data[:, 0]
E = data[:, 1]
M = data[:, 5]
chi = data[:, 6]
CV = data[:, 7]
RunTime = data[:, 8]

plot(T, E,'r',label=r'$40\times40$')
plot(T, data1[:,1],'b',label=r'$60\times60$')
plot(T, data2[:,1],'g',label=r'$80\times80$')
plot(T, data3[:,1],'k',label=r'$100\times100$')
legend(fontsize = l_size)
xlabel(r" Temperature [$kT/J$]", size=f_size)
ylabel("Energy", size=f_size)
xticks(size=f_size-2, rotation=30)
yticks(size=f_size-2, rotation=30)
tight_layout()
show()
plot(T, M,'r',label=r'$40\times40$')
plot(T, data1[:,5],'b',label=r'$60\times60$')
plot(T, data2[:,5],'g',label=r'$80\times80$')
plot(T, data3[:,5],'k',label=r'$100\times100$')
legend(fontsize = l_size)
xlabel(r" Temperature [$kT/J$]", size=f_size)
ylabel("Magnetization", size=f_size)
xticks(size=f_size-2, rotation=30)
yticks(size=f_size-2, rotation=30)
tight_layout()
show()
plot(T, chi,'r',label=r'$40\times40$')
plot(T, data1[:,6],'b',label=r'$60\times60$')
plot(T, data2[:,6],'g',label=r'$80\times80$')
plot(T, data3[:,6],'k',label=r'$100\times100$')
legend(fontsize = l_size, loc='upper right')
xlabel(r" Temperature [$kT/J$]", size=f_size)
ylabel("Susceptibility", size=f_size)
xticks(size=f_size-2, rotation=30)
yticks(size=f_size-2, rotation=30)
tight_layout()
show()
plot(T, CV,'r',label=r'$40\times40$')
plot(T, data1[:,7],'b',label=r'$60\times60$')
plot(T, data2[:,7],'g',label=r'$80\times80$')
plot(T, data3[:,7],'k',label=r'$100\times100$')
legend(fontsize = l_size)
xlabel(r" Temperature [$kT/J$]", size=f_size)
ylabel("Heat capacity", size=f_size)
xticks(size=f_size-2, rotation=30)
yticks(size=f_size-2, rotation=30)
tight_layout()
show()
##############################################################################

# energies = fromfile("energies.bin",dtype=int32)
# indices = int(round(len(energies)*0.9))
# sigma10 = 0.00774658
#
# energies24 = fromfile("energiesT24.bin",dtype=int32)
# indices24 = int(round(len(energies24)*0.9))
# sigma24 = 8.1092
# #
# # x = linspace(-800,-200,1000)
# #
#
# f_size = 28
#
# data = energies24
#
# n,bins,patches=hist(data,45,normed=1, facecolor='green')
# (mu,sigma) = norm.fit(data)
# y = mlab.normpdf(bins,mu,sigma)
# plot(bins,y,'r--',linewidth=2)
# xticks(size=f_size-2, rotation=30)
# yticks(size=f_size-2, rotation=30)
# ylabel(r"Probability", size=f_size)
# xlabel(r"Energy", size=f_size)
# legend(["Gaussian","Energy"], fontsize=f_size-8)
# tight_layout()
# show()
