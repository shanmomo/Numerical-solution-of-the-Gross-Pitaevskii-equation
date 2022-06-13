using BenchmarkTools
using Plots, LaTeXStrings, Measures; gr()
include("TSSP.jl")

xgrid = LinRange(-6,6,800)
tgrid = LinRange(0,1000,50000)
C = 1/24/(xgrid[2]-xgrid[1])^2
function energy(q,Vfunc)
    V = Vfunc.(xgrid,q)
    Aq = C * ( circshift(q,2) - 16*circshift(q,1) + 30*q - 16*circshift(q,-1) + circshift(q,-2) ) + V .* q
    E = real( dot(q,Aq) / dot(q,q) )
end

V(x) = x^2 / 2# + x^4/20 + 4exp(-2x^2)
V0(x,ψ) = V(x)
V1(x,ψ) = V(x) + 5.0 * abs2(ψ)
V2(x,ψ) = V(x) - 2.0 * abs2(ψ)
ψ0 = ComplexF64[ exp(-x^2 /2) / (π)^(1/4) + 0.1randn() for x in xgrid]

ψs0 = iTSSP(ψ0, V0, xgrid, tgrid)
ψs1 = iTSSP(ψ0, V1, xgrid, tgrid)
ψs2 = iTSSP(ψ0, V2, xgrid, tgrid)
E0 = energy( ψs0[:,end-1], V0 )
E1 = energy( ψs1[:,end-1], V1 )
E2 = energy( ψs2[:,end-1], V2 )

p1= plot(size = (800,500), minorgrid = true,
        titlefontsize=16,legendfontsize=12,
        guidefontsize=14,tickfontsize=10,
        legend=:topright,margins=5mm)
    title!("ground-state condensate wavefunction in BEC")
    plot!(xgrid, x->0.05*(V(x).-V(0)), line=(2,0.5), labels="external potential")
    plot!(xgrid, abs2.(ψs0[:,end-1]), line=(4,0.5,:dot), labels="non-interaction")
    plot!(xgrid, abs2.(ψs1[:,end-1]), line=(3,0.5), labels="repulsive interaction")
    plot!(xgrid, abs2.(ψs2[:,end-1]), line=(3,0.5), labels="attractive interaction")
    xaxis!(L"x")
    yaxis!(L"|\psi|^2")
    xlims!(-5,5)
    ylims!(-0.2,0.8)
savefig(p1,"GPE.svg")

p2= plot(size = (800,500), minorgrid = true,
        titlefontsize=16,legendfontsize=12,
        guidefontsize=14,tickfontsize=10,
        legend=:topleft,margins=5mm)
    title!("Degenerate ground state")
    plot!(xgrid, x->0.05*(V(x).-V(0)), line=(2,0.5), labels="external potential")
    plot!(xgrid, abs2.(ψs0[:,end-1]), line=(3,0.5,:dot), labels="non-interaction")
    plot!(xgrid, abs2.(ψs2[:,end-1]), line=(3,0.5), labels="attractive interaction")
    xaxis!(L"x")
    yaxis!(L"|\psi|^2")
    xlims!(-5,5)
    ylims!(-0.2,0.8)
savefig(p2,"Degenerate2.svg")
