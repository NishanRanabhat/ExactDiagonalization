using LinearAlgebra

"""
This file consists of all the bitwise operation functions.

Inputs:

conf -> statevector as integer
i/j -> site indices
N -> system size
"""

function readsite(conf::Int64,i::Int64,N::Int64)
    
    """
    reads the state of i^th site in the state
    """
    
    return (conf&(1<<(N-i)))>>(N-i)
end

function spins_flip(conf::Int64,i::Int64,j::Int64,N::Int64,typeof::String)
    
    """
    
    flips the two sites i and j in a given conf
    
    applies S+ S- + S- S+
    """
    
    if readsite(conf,i,N) != readsite(conf,j,N)
        
        if typeof == "Pauli"
    
            return 2.0, conf ⊻ ((1<<(N-i) ⊻ (1<<(N-j))))
        else
            return 1.0/2.0, conf ⊻ ((1<<(N-i) ⊻ (1<<(N-j))))
        end
    else
        return 0.0, conf
    end
    
end 

function Z(conf::Int64,i::Int64,N::Int64,typeof::String)
    
    """
    applies sZ to i^th site of conf
    """
    
    if typeof == "Pauli"
        return 2.0*readsite(conf,i,N)-1.0, conf
    else
        return readsite(conf,i,N)-1.0/2.0, conf
    end
end

function ZZ(conf::Int64,i::Int64,j::Int64,N::Int64,typeof::String)
    
    """
    applies sZsZ to i^th and j^th sites of conf
    """
    
    si = 2.0*readsite(conf,i,N)-1.0
    sj = 2.0*readsite(conf,j,N)-1.0
    
    if typeof == "Pauli"
        return si*sj, conf
    else
        return si*sj/4.0, conf
    end
end

function X(conf::Int64,i::Int64,N::Int64,typeof::String)

    """
    applies sX to i^th site of conf
    """
    
    if typeof == "Pauli"
        return 1.0, conf ⊻ (1<<(N-i))
    else
        return 1.0/2.0,conf ⊻ (1<<(N-i))
    end
end

function XX(conf::Int64,i::Int64,j::Int64,N::Int64,typeof::String)
    
    """
    applies sXsX to i^th and j^th sites of conf
    """
    
    if typeof == "Pauli"
        return 1.0, X(X(conf,i,N,typeof)[2],j,N,typeof)[2]
    else
        return 1.0/4.0, X(X(conf,i,N,typeof)[2],j,N,typeof)[2]
    end
end


function Y(conf::Int64,i::Int64,N::Int64,typeof::String)

    """
    applies sY to i^th site of conf
    """
    
    if typeof == "Pauli"
        return (2.0*readsite(conf,i,N)-1.0)*im, conf ⊻ (1<<(N-i))
    else
        return (readsite(conf,i,N)-1.0/2.0)*im, conf ⊻ (1<<(N-i))
    end
end

function YY(conf::Int64,i::Int64,j::Int64,N::Int64,typeof::String)
    
    """
    applies sYsY to i^th and j^th sites of conf
    """
    
    si = (2.0*readsite(conf,i,N)-1.0)*im
    sj = (2.0*readsite(conf,j,N)-1.0)*im
    
    if typeof == "Pauli"
        return si*sj, Y(Y(conf,i,N,typeof)[2],j,N,typeof)[2]
    else
        return si*sj/4.0, Y(Y(conf,i,N,typeof)[2],j,N,typeof)[2]
    end
end


function shift_bits(conf::Int64,N::Int64)
    
    """
    translates bits one step to right
    """
    
    return (conf>>1)|(readsite(conf,N,N)<<(N-1))
end

