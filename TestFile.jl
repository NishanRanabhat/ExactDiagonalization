using LinearAlgebra

"""
Note: to use this module you need to download all the files in one folder in your local device
"""

#include the path to the module
include("/home/nishan/CodeDev/ExactDiagonalization/ExactDiag.jl")

#call module
using .ExactDiag

"""
Now you can access all the files and features within this module.
Below is a specific example.  
"""

#Define Parameters
J1=1.0;
J2=1.0
N = 8 #system size

#Define the hilbert space based on choosen parameters
Basis = ExactDiag.hilbert_space(N=N);

#Define index list as needed, here we have taken nearest neighbor interaction with PBC
index_list_one_site = [i for i in 1:N]
index_list_NN = [(i,mod1(i+1,N)) for i in 1:N-1];
index_list_NNN = [(i,mod1(i+2,N)) for i in 1:N-2];

#Define two site operator objects 
OperatorZZ_NN = ExactDiag.Operator_two_site(operator_name = "ZZ",operator_type="Pauli",indices=index_list_NN); 
OperatorXX_NN = ExactDiag.Operator_two_site(operator_name = "XX",operator_type="Pauli",indices=index_list_NN);
OperatorYY_NN = ExactDiag.Operator_two_site(operator_name = "YY",operator_type="Pauli",indices=index_list_NN); 

OperatorX_NN = ExactDiag.Operator_one_site(operator_name = "X",operator_type="Pauli",indices=index_list_one_site); 
OperatorZ_NN = ExactDiag.Operator_one_site(operator_name = "Z",operator_type="Pauli",indices=index_list_one_site); 



OperatorZZ_NNN = ExactDiag.Operator_two_site(operator_name = "ZZ",operator_type="Pauli",indices=index_list_NNN); 
OperatorXX_NNN = ExactDiag.Operator_two_site(operator_name = "XX",operator_type="Pauli",indices=index_list_NNN);
OperatorYY_NNN = ExactDiag.Operator_two_site(operator_name = "YY",operator_type="Pauli",indices=index_list_NNN); 

HamZZ_NN = ExactDiag.Two_site_matrix(OperatorZZ_NN,Basis,N);
HamXX_NN = ExactDiag.Two_site_matrix(OperatorXX_NN,Basis,N);
HamYY_NN = ExactDiag.Two_site_matrix(OperatorYY_NN,Basis,N);

HamX_NN = ExactDiag.One_site_matrix(OperatorX_NN,Basis,N)
HamZ_NN = ExactDiag.One_site_matrix(OperatorZ_NN,Basis,N)


HamZZ_NNN = ExactDiag.Two_site_matrix(OperatorZZ_NNN,Basis,N);
HamXX_NNN = ExactDiag.Two_site_matrix(OperatorXX_NNN,Basis,N);
HamYY_NNN = ExactDiag.Two_site_matrix(OperatorYY_NNN,Basis,N);


#Ham = -1.0*HamXX_NN - 0.35* HamZ_NN 
Ham = J1*(HamXX_NN+HamYY_NN+HamZZ_NN)+J2*(HamXX_NNN+HamYY_NNN+HamZZ_NNN)

#Perform diagonalization
vals = eigvals(Ham);
vecs = eigvecs(Ham);

println(vals)

#the problem might be here, specially in the initial state issue
function coefficient_X(vals,vecs,t,N)

    return vecs * (exp.(-im*vals*t) .* [sum(conj.(vecs[:,j]))/sqrt(2^N) for j in 1:2^N])
end

function coefficient_Z(vals,vecs,t,N)

    return vecs * (exp.(-im*vals*t) .* [vecs[1,j] for j in 1:2^N])
end

function one_site_expectation(Basis,coefficients,site,N)

    val = 0.0

    for conf in Basis.configs

        coeff, new_conf = ExactDiag.X(conf,site,N,"Pauli")

        val += coeff*conj(coefficients[new_conf+1])*coefficients[conf+1]
    end

    return val
end

function time_list(from,to,dt)
    
    t_list = Array{Any,1}(undef,Int((to-from)/dt))
    
    for t in 1:Int((to-from)/dt)
        
        t_list[t] = from + (t-1)*dt
    end
    return t_list
end

#for t in [0.0,0.5,1.0,1.5,2.0]
#println(coefficient_Z(vals,vecs,t,N))
#end

#val = one_site_expectation(Basis,coefficient_X(vals,vecs,10.0,N),2,N)
#println(val)

#for t in time_list(0.0,12.0,0.1)
#val = one_site_expectation(Basis,coefficient_Z(vals,vecs,t,N),6,N)
#println(val)
#end


