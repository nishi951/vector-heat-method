using PlyIO
using LinearAlgebra
using SparseArrays
using ProgressBars
using Debugger

include("ply_utils.jl")

"""
Object containing info necessary to compute the vector heat method
"""
struct VectorHeatMethod
    V::Array{Float64,2}                             # Vertex coordinates
    F::Array{Array{Int64,1},1}                      # Faces
    E::SparseMatrixCSC{Float64, Int64}              # Edge lengths
    A::SparseVector{Float64, Int64}                 # Face areas
    θ::Dict{Int64, SparseMatrixCSC{Float64, Int64}} # Interior angles
    M::SparseMatrixCSC{Float64, Int64}              # Mass associated to vertices
    φ::Dict{Int64, Dict{Int64, Float64}}            # Ordered interior angles
    # Δ::SparseMatrixCSC{Float64, Int64}              # Scalar Laplacian
    Δ::Symmetric{Float64, SparseMatrixCSC{Float64, Int64}}
    # Δ∇::SparseMatrixCSC{ComplexF64, Int64}          # Connection Laplacian
    Δ∇::Hermitian{ComplexF64, SparseMatrixCSC{ComplexF64, Int64}}
    int::Array{Bool, 1}                             # Indicator for interior mesh points
end

function heron(a::Float64, b::Float64, c::Float64)
    s = (a + b + c)/2
    return sqrt(s*(s - a)*(s - b)*(s - c))
end
lawofcosines(a, b, c) = acos((a^2 + b^2 - c^2)/(2*a*b))

function getInteriorVerticesTangentSpaces(V, F, θ)
    N = size(V, 1)
    Ord = Dict(i => Dict() for i in 1:N) # Orderings on adjacent vertices
    Adj = Dict(i => Set{Int64}() for i in 1:N) # Adjacency map
    println("Computing adjacency maps...")
    for (i, j, k) in ProgressBar(F)
        Ord[i][j] = k
        Ord[j][k] = i
        Ord[k][i] = j
        push!(Adj[i], j), push!(Adj[i], k)
        push!(Adj[j], k), push!(Adj[j], i)
        push!(Adj[k], i), push!(Adj[k], j)
    end
    int = zeros(Bool, N)
    for i in 1:N
        int[i] = length(Ord[i]) > 0 ? checkCycle(Ord[i], first(keys(Ord[i]))) : false
    end

    # Holds normalized angles to adjacent vertices
    φ = Dict(i => Dict() for i in 1:N)
    println("Computing tangent space angles...")
    for i in ProgressBar(1:N)
        if !int[i]
            # Boundary points: set 0 for all adjacent vertices
            adj = unique(Adj[i])
            for j in adj
                φ[i][j] = 0
            end
            continue
        end
        adj = []
        norm_angles = []
        curr = first(keys(Ord[i]))
        total_angle = 0   # Expressed as an angle in radians
        while true
            push!(adj, curr)
            push!(norm_angles, total_angle)
            next = Ord[i][curr]
            total_angle += θ[i][curr, next]
            if next == adj[1] # Loop complete
                # Normalize angle
                norm_angles *= 2*π/total_angle
                for (j, angle) in zip(adj, norm_angles)
                    φ[i][j] = angle
                end
                break
            end
            curr = next
        end
    end
    return int, φ
end

"""
Compute relative angle needed to parallel transport a vector,
expressed as a complex exponential.
"""
r(φ, i, j) = exp(im*(φ[j][i] + π - φ[i][j]))
function computeLaplacians(V, F, int, φ, θ)
    N = size(V, 1)
    Δ = spzeros(N, N)
    Δ∇ = spzeros(Complex{Float64}, N, N)
    for (i, j, k) in ProgressBar(F)
        a = cot(θ[i][j, k])
        b = cot(θ[j][k, i])
        c = cot(θ[k][i, j])
        rij = int[i] * r(φ, i, j)
        rjk = int[j] * r(φ, j, k)
        rki = int[k] * r(φ, k, i)
        rji = conj(rij)
        rkj = conj(rjk)
        rik = conj(rki)
        Δ[[i,j,k],[i,j,k]] += -1/2 * [b+c  -c  -b ;
                                      -c   c+a -a ;
                                      -b   -a  a+b]
        # Δ∇[[i,j,k],[i,j,k]] += -1/2 * [b+c    -c*rij -b*rik;
        #                                -c*rji c+a    -a*rjk;
        #                                -b*rki -a*rkj a+b   ]
        Δ∇[[i,j,k],[i,j,k]] += -1/2 * [b+c -c*rji -b*rki;
                                       -c*rij c+a -a*rkj;
                                       -b*rik -a*rjk a+b]
    end
    return Symmetric(Δ), Hermitian(Δ∇)
end

function preprocess(V, F)
    # Blender duplicates vertices for some reason
    println("Removing duplicate vertices...")
    V, F = removeDuplicateVertices(V, F)
    # Ensure consistent face orientation
    println("Checking consistent mesh orientation")
    if !has_consistent_orientation(V, F)
        throw(DomainError("Ply lacks consistent orientation."))
    end
    N = size(V, 1) # Number of vertices
    A = spzeros(size(F, 1)) # Face areas
    θ = Dict(v => spzeros(N, N) for v in 1:N)
    M = zeros(N)
    E = spzeros(N, N)
    println("Computing angles and areas...")
    for (f, (i, j, k)) in ProgressBar(enumerate(F))
        vi, vj, vk = V[i, :], V[j, :], V[k, :]
        ij = norm(vi - vj)
        jk = norm(vj - vk)
        ki = norm(vk - vi)
        E[i,j] = ij
        E[j,i] = ij
        E[j,k] = jk
        E[k,j] = jk
        E[k,i] = ki
        E[i,k] = ki
        # Area of face
        A[f] = heron(ij, jk, ki)
        # Associate 1/3 of the mass to each vertex
        M[i] += 1/3 * A[f]
        M[j] += 1/3 * A[f]
        M[k] += 1/3 * A[f]
        # Interior angles of face (first index is center of angle)
        θ[i][j,k] = lawofcosines(ij, ki, jk)
        θ[i][k,j] = lawofcosines(ij, ki, jk)
        θ[j][i,k] = lawofcosines(ij, jk, ki)
        θ[j][k,i] = lawofcosines(ij, jk, ki)
        θ[k][i,j] = lawofcosines(ki, jk, ij)
        θ[k][j,i] = lawofcosines(ki, jk, ij)
    end
    println("Computing normalized angles for interior vertices...")
    int, φ = getInteriorVerticesTangentSpaces(V, F, θ)
    println("Computing laplacians...")
    Δ, Δ∇ = computeLaplacians(V, F, int, φ, θ)
    M = spdiagm(0 => M)
    return VectorHeatMethod(V, F, E, A, θ, M, φ, Δ, Δ∇, int)
end
