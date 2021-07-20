using Base: Float64
using DelimitedFiles
using SpecialFunctions

include("bound_matrix.jl")

Nexp = 6

Mr = 1 + 2^Nexp
Nz = Mr
MN = Mr * Nz
Mm1 = Mr - 1
Nm1 = Nz - 1
llp = 2 * (Mr + Nz) - 8

Rmin, Rmax, Zmin, Zmax = readdlm("solver.dat")

dr = (Rmax - Rmin) / Mm1
dz = (Zmax - Zmin) / Nm1
alpha = (Rmax - Rmin) / (Rmax + Rmin)
sh = dz / dr
ss = sh^2

println("Rmin = ", Rmin, " Rmax = ", Rmax, " dr = ", dr)
println("Zmin = ", Zmin, " Zmax = ", Zmax, " dz = ", dz)
println("alpha = ", alpha, " sh = ", sh, " ss = ", ss)

ip = Array{Int}(undef, llp)
jp = Array{Int}(undef, llp)
aux = Array{Float64}(undef, llp, llp)

bound_matrix!()

println(aux[22, 32])