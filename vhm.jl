using PlyIO
using PyPlot
# using Plots
using LinearAlgebra
using Statistics

include("vhm_utils.jl")
include("make_flatmesh.jl")
include("visualize_vhm.jl")
# data = "small_bunny"
# data = "bunny_dec_001"
data = "sphere"
ply = load_ply("meshes/$data.ply")
V, F = extractVF(ply)
# V, F = make_flatmesh(10, 10, 0.1)
vhm = preprocess(V, F)

# Specify initial conditions
# Y0 ∈ C^|V|  # Angle of initial vector(s)
# u0 ∈ R^|V|  # Magnitude of initial vector(s)
# ϕ0 ∈ R^|V|  # Indicator for location of initial vector(s)
N = size(vhm.V, 1)
# u0 = zeros(N)
# ϕ0 = zeros(N)

dir = π/3
Y0 = zeros(Complex{Float64}, N)
# Y0[45] = exp(im*dir)
Y0[20] = exp(im*dir)
# Y0[164] = exp(im*0)
Y0[303] = exp(im*(π/4))
# Y0[54] = exp(im*0)
u0 = abs.(Y0)
ϕ0 = Float64.(u0 .!= 0)

# Only work with interior vertices
# M = Diagonal(vhm.M[vhm.int, vhm.int])
# Δ∇ = Hermitian(vhm.Δ∇[vhm.int, vhm.int])
# Δ = Symmetric(vhm.Δ[vhm.int, vhm.int])
M = Diagonal(vhm.M)
L∇ = Hermitian(vhm.Δ∇)
L = Symmetric(vhm.Δ)

Y = zeros(ComplexF64, N)
u = zeros(N)
ϕ = zeros(N)
h = sum(vhm.E)/nnz(vhm.E)
t = h^2  # t = (mean edge length)^2 is a good heuristic

println("Solving system...")
A1 = M - t*L∇
A2 = M - t*L
# Y[vhm.int] = A1 \ Y0[vhm.int]
# u[vhm.int] = A2 \ u0[vhm.int]
# ϕ[vhm.int] = A2 \ ϕ0[vhm.int]
Y = A1 \ Y0
u = A2 \ u0
ϕ = A2 \ ϕ0

X = Y ./ (abs.(Y)) .* (u ./ ϕ)
# X[vhm.int] = Y[vhm.int]./abs.(Y[vhm.int]).*u[vhm.int]./ϕ[vhm.int]
# Set NaN entries to 0
X[isnan.(X)] .= 0
fig = figure()
ax = fig[:gca](projection="3d")
x, y, z = vhm.V[:, 1], vhm.V[:, 2], vhm.V[:, 3]
scatter3D(x, y, z)
plot_edges(ax, vhm.V, vhm.F; color=:lightgrey)
# Initial vector field
Ytrinsic = getExtrinsicVectors(vhm, Y0)
ax = plot_arrow_field(ax, vhm.V, Ytrinsic, 3*h; label_vertices=false, linewidth=4, color=:red)
Xtrinsic = getExtrinsicVectors(vhm, X)
ax = plot_arrow_field(ax, vhm.V, Xtrinsic, h/2; label_vertices=false, linewidth=2)
ax[:axis]("off")
