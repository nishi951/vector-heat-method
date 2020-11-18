using PyPlot

"""
Intrinsic tangent vector x, expressed as a complex exponential r*e^iθ
Tangent bundle φ
Vertex indexed i
Vertex array V, V[i, :] is coordinates of vertex i
"""
function getExtrinsicTangentVector(φ, x, i, V)
    if x == 0
        return zeros(3)
    end
    vi = V[i, :]
    adj = φ[i]  # Maps indices of vertices to angles
    adj = sort(collect(adj), by=x->x[2])
    # n_adj = length(adj)
    # angles = values(adj)
    # vertices = keys(adj)
    # perm = sortperm(values(adj))
    # angles = angles[perm] # Get both into sorted order
    # vertices = vertices[perm]
    # angles = [angles..., 2*π]
    # vertices = [vertices..., vertices[1]] # Add wrapped vertex
    push!(adj, adj[1][1] => 2*π)
    θ = (angle(x) + 2*π) % (2*π)  # Enforce [0, 2π) range
    # for ((j, k), (θj, θk)) in zip(zip(vertices[1:end-1], angles[1:end-1]),
    #                          zip(vertices[2:end], angles[2:end]))
    for ((j, θj), (k, θk)) in zip(adj[1:end-1], adj[2:end])
        if θj <= θ && θ < θk
            λ = (θ - θj)/(θk - θj)# Interpolation parameter
            vj = V[j, :]
            vk = V[k, :]
            vλ = (1 - λ)*vj + λ*vk
            return (vλ - vi)/norm(vλ - vi) * abs(x)
        end
    end
    throw(DomainError("Angle $θ does not lie in [0, 2π]."))
end
function getExtrinsicVectors(vhm, X)
    Xtrinsic = zeros(N, 3)
    println("Computing extrinsic tangent vectors (for plotting)...")
    for (i, x) in enumerate(X)
        Xtrinsic[i, :] = getExtrinsicTangentVector(vhm.φ, x, i, vhm.V)
    end
    return Xtrinsic
end

function plot_arrow_field(ax, V, Xtrinsic, k; label_vertices=false, quiver_kwargs...)
    # fig = figure()
    # ax = fig[:gca](projection="3d")
    x, y, z = V[:, 1], V[:, 2], V[:, 3]
    Xtrinsic_plot = Xtrinsic .* k .* sqrt.(sum(Xtrinsic .^ 2, dims=2))
    u, v, w = Xtrinsic_plot[:, 1], Xtrinsic_plot[:, 2], Xtrinsic_plot[:, 3]
    # scatter3D(x, y, z)
    if label_vertices
        for i in 1:size(x, 1)
            ax[:text](x[i], y[i], z[i], i)
        end
    end
    ax[:quiver](x, y, z, u, v, w; quiver_kwargs...)
    return ax
end

function plot_edges(ax, V, F; plot_kwargs...)
    for (i, j, k) in F
        vi, vj, vk = V[i, :], V[j, :], V[k, :]
        for (r, s) in [(vi, vj), (vi, vk), (vj, vk)]
            ax[:plot]([r[1], s[1]], [r[2], s[2]], [r[3], s[3]]; plot_kwargs...)
        end
    end
    return ax
end
