using LinearAlgebra
using Plots
using PyCall

# Import NumPy via PyCall
np = pyimport("numpy")


# Define spin operators for a single spin
Sx = [0 1; 1 0]
Sy = [0 -im; im 0]
Sz = [1 0; 0 -1]
I = [1 0; 0 1]# Identity matrix


# Utility: Tensor product of multiple matrices
function kron_all_State(st::Any)
    result = st[1]
    for s in st[2:end]
        result = kron(result, s)
    end
    return result
end

function Opt_list(N,j,dir::String)

    List = Array{Any,1}(undef,N)

    if dir == "X"
        unit = [0 1;1 0]
    elseif dir == "Y"
        unit = [0 -im;im 0]
    elseif dir == "Z"
        unit = [1  0;0 -1]
    end

    I = [1 0;0 1]

    for i in 1:N
        if i == j
            List[i] = unit
        else
            List[i] = I
        end
    end

    return List
end

function single_site_Opt(List)

    result = List[1]
    for l in List[2:end]
        result = kron(result, l)
    end
    return result
end

# Utility: Tensor product of multiple matrices
function kron_all_Opt(ops::Any)
    result = ops[1]
    for op in ops[2:end]
        result = kron(result, op)
    end
    return result
end

# Build J1-J2 Hamiltonian
function J1_J2_Hamiltonian(J1, J2, N)
    H = zeros(ComplexF64, 2^N, 2^N)  # Initialize Hamiltonian
    # NN interactions
    for i in 1:(N-1)
        for op in [Sx, Sy, Sz]
            ops = Vector{Matrix{ComplexF64}}()  # Explicitly define as vector of ComplexF64 matrices
            append!(ops, [I for _ in 1:(i-1)])  # Identity matrices before the interaction
            append!(ops, [op, op])             # Add two interacting operators
            append!(ops, [I for _ in (i+2):N])  # Identity matrices after the interaction
            H += J1 * kron_all_Opt(ops)
        end
    end
    # NNN interactions
    for i in 1:(N-2)
        for op in [Sx, Sy, Sz]
            ops = Vector{Matrix{ComplexF64}}()  # Explicitly define as vector of ComplexF64 matrices
            append!(ops, [I for _ in 1:(i-1)])  # Identity matrices before the interaction
            append!(ops, [op, I, op])          # Add two interacting operators with one identity in between
            append!(ops, [I for _ in (i+3):N])  # Identity matrices after the interaction
            H += J2 * kron_all_Opt(ops)
        end
    end
    return H
end

function time_list(from,to,dt)
    
    t_list = Array{Any,1}(undef,Int((to-from)/dt))
    
    for t in 1:Int((to-from)/dt)
        
        t_list[t] = from + (t-1)*dt
    end
    return t_list
end

function fully_pol(N,dir::String,ploarization::String)

    state_list = Array{Any,1}(undef,N)

    if dir == "X"
        if ploarization == "up"
            unit = 1/sqrt(2) * [1; 1]
        elseif ploarization == "down"
            unit = 1/sqrt(2) * [1; -1]
        end
    elseif dir == "Y"
        if ploarization == "up"
            unit = 1/sqrt(2) * [1; im]
        elseif ploarization == "down"
            unit = 1/sqrt(2) * [1; -im]
        end
    elseif dir == "Z"
        if ploarization == "up"
            unit = [1; 0]
        elseif ploarization == "down"
            unit = [0; 1]
        end
    end

    for n in 1:N
        state_list[n] = unit
    end

    return kron_all_State(state_list)
end

function neel(N,dir::String)

    state_list = Array{Any,1}(undef,N)

    if dir == "X"
        unit_up = 1/sqrt(2) * [1; 1]
        unit_down = 1/sqrt(2) * [1; -1]
    elseif dir == "Y"
        unit_up = 1/sqrt(2) * [1; im]
        unit_down = 1/sqrt(2) * [1; -im]
    elseif dir == "Z"
        unit_up = [1; 0]
        unit_down = [0; 1]
    end

    for n in 1:N
        if mod(n, 2) == 0
            state_list[n] = unit_up
        else
            state_list[n] = unit_down
        end
    end

    return kron_all_State(state_list)
end

function evol_expect(J1,J2,N,j,init_state,Operator,Ham)

    count = 0
    avg_list = Array{Any,1}(undef,2000)
    for t in time_list(0.02,40.02,0.02)
        count += 1
        state = exp(-im*t*Ham)*init_state
        avg = state'*Operator*state
        avg_list[count] = avg
        println(t,",",avg)
    end
    return avg_list
end

J1 = 1.0
J2 = 1.0
N = 9
dir = "X"
j = 5
#Ham = J1_J2_Hamiltonian(J1,J2,N)
#vals = eigvals(Ham);
#vecs = eigvecs(Ham);

init_state = neel(N,dir)
List = Opt_list(N,j,"Z")
Operator = single_site_Opt(List)
Ham = J1_J2_Hamiltonian(J1, J2, N)

avg_list = evol_expect(J1,J2,N,j,init_state,Operator,Ham)

avg_list = np.array(avg_list)

np.save("my_array_state_X_opt_Z.npy", avg_list)




