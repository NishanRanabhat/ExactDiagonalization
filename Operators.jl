using LinearAlgebra
using Parameters

"""
The important idea is to create the Operator_one/two_site object which acts like an operator as we would use in text book. 

fields :

operator_name -> locates the corresponding operator function.
indices -> set of indices where the operator acts. 
           This way we can create any kind of interaction among two sites by manually creating interaction between two spins of arbritrary kind.
           
Furthermore, we can also add another field in this struct (or just add another element in the indices field) that will give the coupling 
constant corresponding to the given operator and indices. This way we can create any kind of spin hamiltonian upto two spin interaction.
"""

@with_kw struct Operator_two_site
    
    operator_name::AbstractString
    operator_type::AbstractString = "SpinHalf"
    indices::Vector{Tuple{Int64, Int64}}
end 

@with_kw struct Operator_one_site
    
    operator_name::AbstractString
     operator_type::AbstractString = "SpinHalf"
    indices::Vector{Int64}
end   

function Two_site_matrix(operator::Operator_two_site,hilbert_space::Union{basis1,basis2},N::Int64)
    
    Ham = zeros(length(hilbert_space.configs),length(hilbert_space.configs))
        
    count = 0
    
    oprt = eval(Symbol(operator.operator_name))
    oprt_typ = operator.operator_type
        
    if nfields(hilbert_space) == 1
        
        for conf in hilbert_space.configs
        
            count += 1
        
            for indx in operator.indices
                            
                coeff, new_conf = oprt(conf,indx[1],indx[2],N,oprt_typ)
                            
                Ham[bisection_checkstate(hilbert_space,new_conf),count] += coeff
            end
        end
        
        return Ham
        
    elseif nfields(hilbert_space) == 3
        
        for conf in hilbert_space.configs
        
            count += 1
        
            for indx in operator.indices
            
                coeff, new_conf = oprt(conf,indx[1],indx[2],N,oprt_typ)
                
                represent, l = representative(new_conf,N)
                
                b = bisection_checkstate(hilbert_space,represent)
                
                if b >= 0
                
                    Ham[b,count] += coeff_trans_inv(hilbert_space.periodicity[count],hilbert_space.periodicity[b],coeff,hilbert_space.k,l,N)
                end                
            end
        end
        
        return Ham
    end
end

function One_site_matrix(operator::Operator_one_site,hilbert_space::Union{basis1,basis2},N::Int64)
    
    Ham = zeros(length(hilbert_space.configs),length(hilbert_space.configs))
        
    count = 0
    
    oprt = eval(Symbol(operator.operator_name))
    oprt_typ = operator.operator_type
    
    if nfields(hilbert_space) == 1
        
        for conf in hilbert_space.configs
        
            count += 1
        
            for indx in operator.indices
            
                coeff, new_conf = oprt(conf,indx,N,oprt_typ)
            
                Ham[bisection_checkstate(hilbert_space,new_conf),count] += coeff
            end
        end
        
        return Ham
        
    elseif nfields(hilbert_space) == 3
        
        for conf in hilbert_space.configs
        
            count += 1
        
            for indx in operator.indices
            
                coeff, new_conf = oprt(conf,indx,N,oprt_typ)
                
                represent, l = representative(new_conf,N)
                
                b = bisection_checkstate(hilbert_space,represent)
                
                if b >= 0
                
                    Ham[b,count] += coeff_trans_inv(hilbert_space.periodicity[count],hilbert_space.periodicity[b],coeff,hilbert_space.k,l,N)
                end                
            end
        end
        
        return Ham 
    end
end