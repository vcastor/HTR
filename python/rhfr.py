#!/opt/homebrew/bin/python3
#-----------------------------------------------------------------------
#***********************************************************************
#  =This program do the Rothaan method for Hartree-Fock calculation=
#***********************************************************************
# The same as all versions... BUT IN python 3 🥳
#***********************************************************************
#       with love,
#               Victoria Castor 2022
#-----------------------------------------------------------------------

import numpy as np
import sys
import math
import scipy as sp
from scipy.linalg import eigh

pi = 3.1415926535
#-----------------------------------------------------------------------
#                             functions
def f_cou(q,Z,dis):
    ecou = q*Z/dis
    return ecou

def f_distance2(r1,r2):
    d1 = r1[0] - r2[0]
    d2 = r1[1] - r2[1]
    d3 = r1[2] - r2[2]
    dis2 = d1*d1 + d2*d2 + d3*d3
    return dis2

def f_cero(x):
    if (x > 0.0000000001):
        y = 0.5*np.sqrt(pi/x)*math.erf(np.sqrt(x))
    else:
        y = 1
    return y

#-----------------------------------------------------------------------
#                           READ the input file
inp = sys.argv[1]
inpfile = open(inp, 'r')
datainp = inpfile.readlines()
basis = datainp[0]
q = datainp[1]
q = int(q)
unit = datainp[2]
na = datainp[3]
na = int(na)
xyz = np.zeros((3,na))
atom = []
for i in range(na):
        linen = datainp[i+5]
        linen = linen.split()
        atom.append(linen[0])
        xyz[0][i] = linen[1]
        xyz[1][i] = linen[2]
        xyz[2][i] = linen[3]
inpfile.close()

if unit.lower() == 'angstrom\n':
        xyz = np.multiply(xyz, 1.8897259886)
#elif unit.lower() == 'bohr\n':
#        xyz = xyz

#-----------------------------------------------------------------------
#                                   At basis
#--

general_data = "../fortran/basis/" + basis.rstrip() + "/data_of_basis.dat"

data_basis = open(general_data, 'r')
max_val = data_basis.readlines()
max_bf  = int(max_val[0])
max_pri = int(max_val[1])
data_basis.close()

max_zeta = na*max_bf

file_basis_set = []
for i in range(na):
        file_basis_set.append("../fortran/basis/" + basis.rstrip().lower() + "/" + atom[i].rstrip() + ".dat")

a_num    = np.zeros(na, dtype=int)
n_bf_pa  = np.zeros(na, dtype=int)
n_pri_bf = np.zeros(max_zeta, dtype=int)
zeta     = np.zeros((max_zeta,max_pri))
d        = np.zeros((max_zeta,max_pri))

bf = 0
m  = 0
for i in range(na):
        l = 0
        data_reading = open(file_basis_set[i], 'r')
        reading = data_reading.readlines()
        line_bs = reading[0]
        line_bs = line_bs.split()
        a_num[i]   = int(line_bs[0])
        n_bf_pa[i] = int(line_bs[1])
        l += 1
        bf += n_bf_pa[i]
        for j in range(n_bf_pa[i]):
                n_pri_bf[j+m] = int(reading[l])
                l += 1
                for k in range(n_pri_bf[j+m]):
                        line_zd    = reading[l]
                        line_zd    = line_zd.split()
                        zeta[j+m][k] = float(line_zd[0])
                        d[j+m][k]    = float(line_zd[1])
                        l += 1
        m += n_bf_pa[i]
        data_reading.close()

#-----------------------------------------------------------------------
#                               Distances

xyz_bf  = np.zeros((3,bf))
dis_a   = np.zeros((na,na))
dis2_a  = np.zeros((na,na))
dis_bf  = np.zeros((bf,bf))
dis2_bf = np.zeros((bf,bf))

for i in range(na):
        for j in range(i, na):
                dis2_a[i][j] = f_distance2(xyz[:,i], xyz[:,j])
                dis2_a[j][i] = dis2_a[i,j]
                dis_a[i][j]  = np.sqrt(dis2_a[i][j])
                dis_a[j][i]  = dis_a[i][j]

k = 0
for i in range(na):
        for j in range(n_bf_pa[i]):
                xyz_bf[:,k] = xyz[:,i]
                k += 1

n = 0
m = 0
for i in range(na):
    for j in range(na):
        for k in range(n_bf_pa[i]):
            for l in range(n_bf_pa[j]):
                dis2_bf[n][m] = dis2_a[j][i]
                dis_bf[n][m]  = dis_a[j][i]
                n += 1
            n -= n_bf_pa[j]
            m += 1
        n += n_bf_pa[j]
        m -= n_bf_pa[i]
    m += n_bf_pa[i]
    n = 0

#-----------------------------------------------------------------------
#                            Born-Oppenheimer
#--

Enn = 0.0
Enn_mat = np.zeros((na,na))

for i in range(na-1):
        for j in range(i+1, na):
                Enn_mat[i][j] = f_cou(a_num[i], a_num[j], dis_a[i][j])
                Enn_mat[j][i] = Enn_mat[i][j]
                Enn += Enn_mat[i][j]

#-----------------------------------------------------------------------
#                         Electrons and robitals
#--
n_ele = 0
for i in range(na):
        n_ele += a_num[i]

n_ele -= q
n_occ = n_ele/2
n_occ = int(n_occ)

#-----------------------------------------------------------------------
#                                Integrals
#--

zeta_escalar = np.zeros((max_zeta,max_pri,max_zeta,max_pri))
xi_escalar = np.zeros((max_zeta,max_pri,max_zeta,max_pri))

for i in range(bf):
    for j in range(i,bf):
        for k in range(n_pri_bf[i]):
            for l in range(n_pri_bf[j]):
                zeta_escalar[i][k][j][l] = zeta[i][k] + zeta[j][l]
                xi_escalar[i][k][j][l]   = zeta[i][k] * zeta[j][l] / zeta_escalar[i][k][j][l]

kfact = 2.0*pi**(5.0/2.0)

rp = np.zeros(3)
rq = np.zeros(3)
S = np.zeros((bf,bf))
T = np.zeros((bf,bf))
V = np.zeros((bf,bf))
Two_ele = np.zeros((bf,bf,bf,bf))

for i in range(bf):
    for j in range(i,bf):
        for k in range(n_pri_bf[i]):
            for l in range(n_pri_bf[j]):
                rp = (zeta[i][k]*xyz_bf[:,i] + zeta[j][l]*xyz_bf[:,j])/zeta_escalar[i][k][j][l]
                overl = np.exp(-xi_escalar[i][k][j][l]*dis2_bf[i][j]) * (pi/zeta_escalar[i][k][j][l])**(3.0/2.0)
                S[i][j] += d[i][k]*d[j][l]*overl
                kine  = xi_escalar[i][k][j][l]*(3.0-(2.0*xi_escalar[i][k][j][l]*dis2_bf[i][j]))*overl
                T[i][j] += d[i][k]*d[j][l]*kine
                NucEle = 0.0
                for m in range(na):
                    rpm2 = f_distance2(rp, xyz[:,m])
                    NucEle += 2.0*a_num[m]*np.sqrt(zeta_escalar[i][k][j][l]/pi)*overl*f_cero(zeta_escalar[i][k][j][l]*rpm2)
                V[i][j] -= d[i][k]*d[j][l]*NucEle
                Kkl = (1.0/zeta_escalar[i][k][j][l])*np.exp(-xi_escalar[i][k][j][l]*dis2_bf[i][j])
                for m in range(bf):
                    for n in range(m,bf):
                        if (int(j*(j+1)/2 +i) >= int(n*(n+1)/2 +m)):
                            for o in range(n_pri_bf[m]):
                                for p in range(n_pri_bf[n]):
                                    rq = (zeta[m][o]*xyz_bf[:,m] + zeta[n][p]*xyz_bf[:,n])/zeta_escalar[m][o][n][p]
                                    Kmn = (1.0/zeta_escalar[m][o][n][p])*np.exp(-xi_escalar[m][o][n][p]*dis2_bf[m][n])
                                    rho = (zeta_escalar[i][k][j][l]*zeta_escalar[m][o][n][p])/(zeta_escalar[i][k][j][l] + zeta_escalar[m][o][n][p])
                                    rpq2 = f_distance2(rp,rq)
                                    phi = (kfact*Kkl*Kmn/np.sqrt(zeta_escalar[i][k][j][l]+zeta_escalar[m][o][n][p]))*f_cero(rho*rpq2)
                                    Two_ele[i][j][m][n] += d[i][k]*d[j][l]*d[m][o]*d[n][p]*phi

for i in range(bf):
    for j in range(i,bf):
        if (i != j):
            S[j][i] = S[i][j]
            V[j][i] = V[i][j]
            T[j][i] = T[i][j]
        for k in range(bf):
            for l in range(k,bf):
                if (int(j*(j+1)/2 + i) >= int(l*(l+1)/2 +k)):
                    Two_ele[i][j][l][k] = Two_ele[i][j][k][l]
                    Two_ele[j][i][k][l] = Two_ele[i][j][k][l]
                    Two_ele[j][i][l][k] = Two_ele[i][j][k][l]
                    Two_ele[k][l][i][j] = Two_ele[i][j][k][l]
                    Two_ele[k][l][j][i] = Two_ele[i][j][k][l]
                    Two_ele[l][k][i][j] = Two_ele[i][j][k][l]
                    Two_ele[l][k][j][i] = Two_ele[i][j][k][l]
    S[i][i] = 1.0

H = T + V
C = np.zeros((bf,bf))
E_old = -1.0
converged=False
step = 0
while(not converged):

    P = np.zeros((bf,bf))
    for i in range(bf):
        for j in range(bf):
            for k in range(n_occ):
                    P[i][j] += C[i][k]*C[j][k]

    J = np.zeros((bf,bf))
    for i in range(bf):
        for j in range(bf):
            for k in range(bf):
                for l in range(bf):
                    J[i][j] += P[k][l]*Two_ele[i][j][l][k]
                    
    K = np.zeros((bf,bf))
    for i in range(bf):
        for j in range(bf):
            for k in range(bf):
                for l in range(bf):
                    K[i][j] += P[k][l]*Two_ele[i][l][k][j]         

    F = H + 2*J - K
    
    # FC=SCE
    E,C = eigh(F, S, eigvals_only=False)
    
    # Energy
    E_elec = 0.0
    for i in range(bf):
        for j in range(bf):
            E_elec = E_elec + P[j][i]*(H[i][j] + F[i][j])

    
    #print("Electronic energy: ", E_elec)

    step += 1
    #converged?
    if(abs(E_old - E_elec)< 0.0000001 or step > 100):
        converged = True
    else:
        E_old = E_elec
        
print('Electronic Energy', E_elec)
print('NN Energy', Enn)
E_T = E_elec + Enn
print("Total energy: ", E_T)

