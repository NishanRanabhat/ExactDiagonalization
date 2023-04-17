using LinearAlgebra
import Base.@kwdef
using KeywordDispatch

"""
We create two kinds of basis objects;

basis1 -> for systems with no translation invariance symmetry implemented. 

          field :

          configs -> which is the vector of states represented as integer digits.

basis2 -> for systems with translation invariance implemented. 

          fields: 

          configs -> vector of states represented as integer digits, 
          preiodicity -> vector of integer digits representing the periodicity of the states in field configs, 
          k -> real number representing the momentum sector.
"""

struct basis1
    
    configs::Vector{Int64} 
end

struct basis2
    
    configs::Vector{Int64}  
    periodicity::Vector{Int64}
    k :: Real
end 

"""
Different methods instances of the function hilbert_space using multiple dispatch.

Inputs:

N -> system size
Nup -> number of up spins
k -> momentum quantum number
"""

@kwdispatch hilbert_space()

@kwmethod function hilbert_space(;N::Int64)
    
    hil_spc = zeros(Int64,2^N)
    
    count = 0
    
    for i in 0:2^(N)-1
        
        count += 1
        
        hil_spc[count] = i
    end
    
    return basis1(hil_spc)
end


@kwmethod function hilbert_space(;N::Int64,Nup::Int64)
    
    if Nup > N
        error("Nup must be less than or equal to N")
    else
    
        hil_spc = zeros(Int64,binomial(N,Nup))
    
        count = 0
    
        for i in 0:2^(N)-1
                        
            if count_ones(i) == Nup
            
                count += 1
            
                hil_spc[count] = i
            end
        end
        return basis1(hil_spc)
    end
end

@kwmethod function hilbert_space(;N::Int64,k::Real)
        
    hil_spc = zeros(Int64,(2^N,2))
    
    count = 0
            
    for conf in hilbert_space(N=N).configs
                
        r = checkstate(conf,k,N)
        
        if r >= 0
            
            count += 1

            hil_spc[count,1] = conf
            hil_spc[count,2] = r
        end
    end
    
    return basis2(hil_spc[1:count,1],hil_spc[1:count,2],k)
end 

@kwmethod function hilbert_space(;N::Int64,Nup::Int64,k::Real)
    
    hil_spc = zeros(Int64,(2^N,2))
    
    count = 0
            
    for conf in hilbert_space(N=N,Nup=Nup).configs
                
        r = checkstate(conf,k,N)
        
        if r >= 0
            
            count += 1

            hil_spc[count,1] = conf
            hil_spc[count,2] = r
        end
    end
    
    return basis2(hil_spc[1:count,1],hil_spc[1:count,2],k)
end 
