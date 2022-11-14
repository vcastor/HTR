#


# functions

function f_cou()

end

f_distance2


f_cero


# READ input file

# basis set

# distances

xyz_bf  = zeros((3.bf))
dis_a   = zeros((na,na))
dis2_a  = zeros((na,na))
dis_bf  = zeros((bf,bf))
dis2_bf = zeros((na,na))

for i in 1:na
        for j in i:na
                dis2_a[i,j] = f_distance2(xyz[:,i],xyz[:,j])
                dis2_a[j,i] = dis2_a[i,j]
                dis_a[i,j]  = sqrt(dis2_a[i,j])
                dis_a[j,i]  = dis_a[i,j]

k = 0
for i in 1:na
        for j in 1:n_bf_pa[i]
                xyz_bf[:,k] = xyz[:,i]
                k += 1
# Borh-Oppenheimer

# Electrons and orbitals

# Integrals

# SCF algorithm HARTREE FOCK ROOTHAAN

# print the result 
