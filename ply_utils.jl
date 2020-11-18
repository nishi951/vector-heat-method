using PlyIO

function extractVF(ply::Ply)
    V = cat((ply["vertex"][d] for d in ["x", "y", "z"])..., dims=2)
    V = Float64.(V)
    F = Array(ply["face"]["vertex_indices"])
    F = [Int64.(f) .+ 1 for f in F]  # Faces are 0-indexed
    return V, F
end

indexof(elem, arr) = findfirst(i -> i == elem, arr)
function removeDuplicateVertices(V, F)
    newaxis = [CartesianIndex()]
    N = size(V, 1)
    remaining = Set(1:N)
    # synonyms maps old vertex indices to new vertex indices
    synonyms = zeros(Int64, N)
    while length(remaining) > 0
        i = pop!(remaining)
        synonyms[i] = i
        v = V[i, :]
        dupes = Set(findall(all(v[newaxis, :] .== V, dims=2)[:]))
        for d in dupes
            synonyms[d] = i
        end
        setdiff!(remaining, dupes)
    end
    new_verts = unique(synonyms) # New list
    V = V[new_verts, :]
    for (j, f) in enumerate(F)
        F[j] = map(i -> indexof(synonyms[i], new_verts), f)
    end
    return V, F
end

function has_consistent_orientation(V, F)
    N = size(V,1)
    # adjacent vertices should always appear in the opposite order across faces
    # or only once in one direction, over all faces
    E = Dict(v => [] for v in 1:N)
    for (i, j, k) in ProgressBar(F)
        push!(E[i], j)
        push!(E[j], k)
        push!(E[k], i)
    end

    for i in 1:N
        # @assert length(E[i]) == length(unique(E[i])) # No edge added twice in same orientation
        # return false
        if length(E[i]) != length(unique(E[i]))
            return false
        end
    end
    return true
end

function checkCycle(d::Dict, start)
    L = 0 # Length of cycle
    curr = start
    visited = Set()
    while true
        if !(d[curr] in keys(d)) || curr in visited
            return false
        end
        push!(visited, curr)
        curr = d[curr]
        L += 1
        if curr == start && L == length(d)
            return true
        end
    end
end
