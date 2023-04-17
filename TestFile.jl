using LinearAlgebra
using Arpack

"""
Note: to use this module you need to download all the files in one folder in your local device
"""

#include the path to the module
include("/path/to/module/folder/ExactDiag.jl")

#call module
using .ExactDiag

"""
Now you can access all the files and features within this module.
Below is a specific example.  
"""

#Define Parameters

N = 8 #system size
Nup = 4 #number of up spins,

#Define the hilbert space based on choosen parameters
Basis = ExactDiag.hilbert_space(N=N,Nup=Nup);

#Define index list as needed, here we have taken nearest neighbor interaction with PBC
index_list = [(i,mod1(i+1,N)) for i in 1:N];

#Define two site operator objects 
OperatorZZ = ExactDiag.Operator_two_site(operator_name = "ZZ",operator_type="Pauli",indices=index_list); #ZZ interaction
OperatorXY = ExactDiag.Operator_two_site(operator_name = "spins_flip",operator_type="Pauli",indices=index_list); # (XX+YY) interaction

#Build Hamiltonian matrices for given operators and Basis
HamZZ = ExactDiag.Two_site_matrix(OperatorZZ,Basis,N);
HamXY = ExactDiag.Two_site_matrix(OperatorXY,Basis,N);
Ham = HamZZ+HamXY;

#Perform diagonalization
vals = eigvals(Ham);
vecs = eigvecs(Ham);