using LinearAlgebra

"""
This file consists helper functions for the package 
"""

function checkstate(conf::Int64,k::Real,N::Int64)
        
    """
    this function checks if a given conf is a valid representation 
    and if it is returns the periodicity, also checks if the conf 
    is compatible with the given momentum.

    Inputs:

    Conf -> statevector as integer
    k -> momentum quantum number
    N -> system size
    """

    R = -1; t = conf
    
    for i in 1:N
        
        t = shift_bits(t,N)
        
        if t < conf 
            return -1
        elseif t == conf
            if mod(k,N/i) != 0 
                return -1
            else
                return i
            end
        end
    end
end


function bisection_checkstate(hilbert_space,conf::Int64)

    """
    checks the index corresponding to a configuration using bisectional algorithm

    Inputs:

    hilbert_space -> space of search
    conf -> statevector as integer whose index in hilbert_space is to be searched
    """

    bmin = 1; bmax = length(hilbert_space.configs)
    
    count = 0

    while true
        
        count += 1
    
        b = bmin + floor(Int64,(bmax-bmin)/2)
    
        if conf < hilbert_space.configs[b]
        
            bmax = b-1
        
        elseif conf > hilbert_space.configs[b]
        
            bmin = b + 1
        
        else
            return b
        
            break
        end
        
        if bmin > bmax
            
            b = -1
            
            return b
        end
        
        if count > floor(Int64,length(hilbert_space.configs)/2)+1
            
            println("Configuration doesn't exist in this Hilbert space")
            
            break
        end
        
    end
end


function representative(new_conf::Int64,N::Int64)
    

    """
    Takes a configuration and returns its representative i.e. the smallest integer state in its 
    equivalent class and the number of shift operators required to reach that state aka the distance
    between a state and its representative.

    Inputs:

    new_conf -> statevector as integer whose representative is to be determined
    N -> system size
    """

    t = new_conf
    represent = new_conf
    l = 0
    
    for i in 1:N
        
        t = shift_bits(t,N)
                        
        if t < represent
            
            represent = t
            l = i
        end
    end
    
    return represent, l
end 


function coeff_trans_inv(Ra::Int64,Rb::Int64,coeff::Number,k::Real,l::Int64,N::Int64)

    """
    gives coefficient corresponding to translation invariant matrix
    """
    
    return coeff*(Ra/Rb)^(0.5)*exp(2.0*pi*k*l/N)
end


"""
Helper functions to connecting string <-> integer <-> vector representation of a statevector
"""

function integer_to_string(n::Int64,N::Int64)
    
    """
    Converts given integer to bitwise string
    n : integer input
    N : size of spin chain
    """
    
    if n <= 2^(N)-1
    
        return bitstring(n)[64-N+1:64]
        
    else
        println("Error: integer bigger than what present number of bits can reprsent")
    end
end

function string_to_integer(str::String)
    
    """
    converts bitwise string to integer
    """
    
    return parse(Int,str; base=2)
end

function string_to_vector(str::String)
    
    """
    converts bitwise string to vector of 0 and 1
    """
    
    return [parse(Int,str[i]) for i in 1:length(str)]
end

function vector_to_string(state::Vector{Int})
    
    """
    converts vector of 0 and 1 to bitwise string
    """
        
    str = string(state[1])
    
    for i in 2:length(state)
        
        str = string(str,state[i])
    end

    return str
end     

function vector_to_integer(state::Vector{Int})
    
    """
    converts vector of 0 and 1 to integer
    """
    
    return parse(Int,vector_to_string(state),base=2)
end

function integer_to_vector(n::Int64,N::Int64)
    
    """
    converts integer to vector of 0 and 1
    """
    
    return string_to_vector(integer_to_string(n,N))
end 