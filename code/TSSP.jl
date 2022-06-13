## TSSP Method for 1D Non-linear Schrödinger Equation
# ψ0, x should be 1D arrays of the same size.
# V(x,ψ) is the non-linear potential.
using LinearAlgebra
using FFTW

@inline function iV_step!(ψn, V, x, dt, Nx)
    @inbounds for i in 1:Nx
        ψn[i] *= exp(-V(x[i],ψn[i])*dt)
    end
end

@inline function iT_step!(ψn, ω_g, dt, Nx)
    fft!(ψn)
    @inbounds for i in 1:Nx
        ψn[i] *= exp(-ω_g[i]*dt)
    end
    ifft!(ψn)
end

function iTSSP(ψ0, V, x, t)
    dx = x[2] - x[1]
    Nx = length(x); Nt = length(t)
    # k空间的离散点
    k_g = (2π / dx) .* (((-Nx÷2)//Nx) : (1//Nx) : ((Nx÷2-1)//Nx))
    # 计算对应频率并且按照fft的顺序重排
    # ω_g = (k_g.^2/2)
    ω_g = circshift((k_g.^2/2) , -(Nx÷2))
    ψs = zeros(ComplexF64,Nx,Nt) # 分配数组
    ψs[:,1] .= ψ0 # 初值
    for i in 1:Nt-1
        normalize!(@view ψs[:,i])
        ψs[:,i] ./= √dx
        dt = t[i+1]-t[i]
        ψn = @view ψs[:,i+1]
        copy!(ψn,ψs[:,i])
        # Step V-T-V
        iV_step!(ψn, V, x, dt/2, Nx)
        iT_step!(ψn, ω_g, dt, Nx)
        iV_step!(ψn, V, x, dt/2, Nx)
    end
    return ψs
end
