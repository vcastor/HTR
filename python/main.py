#-----------------------------------------------------------------------
#***********************************************************************
#  =This program do the Rothaan method for Hartree-Fock calculation=
#***********************************************************************
# PLEASE DO NOT ATTEMP TO SIMPLIFY THIS CODE.
# KEEP THE SPACE SHUTTLE FLYING.
#***********************************************************************
# We call this style 'space shuttle style'. Space shuttle style is meant
# to ensure that every branch and condition is considered and accounted
# for the same way code is written at NASA for apllications like the --
# space shuttle.
#***********************************************************************
# The example inputs are in the same directory as these code. To run ---
# these code we used python 3
#***********************************************************************
#       with love,
#               VCastor 2021
#-----------------------------------------------------------------------

import psi4
import numpy as np
from scipy.linalg import eigh
import sympy as sp
sp.init_printing()

def fdis2(x1, y1, z1, x2, y2, z2):
        d2 = (x1-x2)**2 + (y1-y2)**2 + (z1-z3)**2
        return d2

def fcoulomb(q, Z, distance):
        coulombnn = q*Z/distance
        return coulombnn

#--
#READ the input file
#--
inpfile = open(inp, 'r')
datainp = inpfile.readlines()
basis = datainp[0]
q = datainp[1]
q = int(q)
unit = datainp[2]
na = datainp[3]
na = int(na)
atom = []
x = []
y = []
z = []
for i in range(na):
        linen = datainp[i+5]
        linen = linen.split()
        atom.append = linen[0]
        x.append    = linen[1]
        y.append    = linen[2]
        z.append    = linen[3]
inpfile.close()

#--
#At basis
#--

#--
#Distances
#--

for i in range(na):
        for j in range(na):
                dis2_a[i,j] = fdis2(x[i], y[i], z[i], x[j], y[j], z[j])
                dis2_a[j,i] = dis2_a[i,j]
                dis_a[i,j]  = np.sqrt(dis2_a[i,j])

#--
#Bron-Oppenheimer
#--

#--
#Electrons and robitals
#--
n_ele = 0
for i in range(na):
        n_ele += a_num[i]

n_ele -= q
n_occ = n_occ/2

#--
#Integrals
#--

S = np.asarray(mints.ao_overlap())
print("----------------Matriz S----------------")
sp.pprint(sp.Matrix(S))

T = np.asarray(mints.ao_kinetic())
V = np.asarray(mints.ao_potential())
H=T+V
print("----------------Matriz H----------------")
sp.pprint(sp.Matrix(H))

I = np.asarray(mints.ao_eri())
#print("Integrales (ij|kl):")
#print(I)

nbf = S.shape[0]
ndocc = wfn.nalpha()
#print(nbf)
#print(ndocc)

C = np.zeros((nbf,nbf))
print(C)

E_old = -1.0
converged=False

while(not converged):
    print("---------------------------------ITERACION---------------------------------")  

    print("------------Matriz C Entrada------------")
    sp.pprint(sp.Matrix(C))    
    
    #Paso 4. Calcular P, J y K
    P = np.zeros((nbf,nbf))
    for i in range(nbf):
        for j in range(nbf):
            for k in range(ndocc):
                    P[i][j] = P[i][j] + C[i][k]*C[j][k]
    print("----------------Matriz P----------------")
    sp.pprint(sp.Matrix(P))
                    
    J = np.zeros((nbf,nbf))
    for i in range(nbf):
        for j in range(nbf):
            for k in range(nbf):
                for l in range(nbf):
                    J[i][j] = J[i][j] + P[k][l]*I[i][j][l][k]
    print("----------------Matriz J----------------")
    sp.pprint(sp.Matrix(J))
                    
    K = np.zeros((nbf,nbf))
    for i in range(nbf):
        for j in range(nbf):
            for k in range(nbf):
                for l in range(nbf):
                    K[i][j] = K[i][j] + P[k][l]*I[i][l][k][j]         
    print("----------------Matriz K----------------")
    sp.pprint(sp.Matrix(K))                    

    #Paso 5. Calcular F = H + 2J - K
    F = H + 2*J - K
    print("----------------Matriz F----------------")
    sp.pprint(sp.Matrix(F))
    
    #Paso 6. Resolver FC=SCE
    E,C = eigh(F, S, eigvals_only=False)
    print("-------------Matriz C Salida-------------")
    sp.pprint(sp.Matrix(C))
    
    #Paso 7. Calcular E=sum_i sum_j P_ji (H_ij+F_ij)
    E_elec = 0.0
    for i in range(nbf):
        for j in range(nbf):
            E_elec = E_elec + P[j][i]*(H[i][j] + F[i][j])

    
    print("Energia Electronica: ",E_elec)

    #Paso 8. ¿$E_i=E_{i-1}$?, Si: acabe. No: volver a paso 4.
    if(abs(E_old - E_elec)< 0.0000001):
        converged = True
    else:
        E_old = E_elec
        
E_nuc = H2.nuclear_repulsion_energy() 
print("Energia nuclear: ", E_nuc)
E_T = E_elec + E_nuc
print("Energia Total: ", E_T)

