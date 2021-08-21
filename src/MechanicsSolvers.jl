@views function StokesSolver!(Vx::Matrix{Float64}, Vy::Matrix{Float64}, Pc::Matrix{Float64}, NumVx::Matrix{Int64}, NumVy::Matrix{Int64}, NumP::Matrix{Int64}, Fx::Matrix{Float64}, Fy::Matrix{Float64}, Fp::Matrix{Float64}, Kuu::SparseMatrixCSC{Float64, Int64}, Kup::SparseMatrixCSC{Float64, Int64}, Kpu::SparseMatrixCSC{Float64, Int64}, etac::Matrix{Float64}, gamma::Float64, solver::Int64)
    fu = zeros(Float64, length(Fx) + length(Fy))
    fp = zeros(Float64, length(Fp))
    fu[1:length(Fx)]     .= Fx[:]
    fu[length(Fx)+1:end] .= Fy[:]
    fp .= Fp[:]
    if solver == 0
        # Slightly compressible pressure block
        cP   = 1.0./gamma.*ones(Float64,size(etac,1), size(etac,2))
        I    = NumP[:]
        J    = NumP[:]
        V    = cP[:]
        PP   = sparse( I, J, V)
        Kmat = [Kuu Kup; Kpu PP]
        F    = [fu; fp]
        dX   = -Kmat\F
        Vx[:,2:end-1] .+= dX[NumVx]
        Vy[2:end-1,:] .+= dX[NumVy]
        Pc            .+= dX[NumP.+maximum(NumVy)]
        Pc            .-= mean(Pc)
    elseif solver == 1
        DecoupledSolver!( Vx,Vy,Pc,NumVx,NumVy,NumP,fu,fp,Kuu,Kup,Kpu,etac,gamma )
    elseif solver == 2 
        @tturbo coef  = gamma*ones(length(etac))#.*etac[:]
        @tturbo Kppi  = spdiagm(coef)
        @tturbo M     = [Kuu Kup; Kpu 0*Kppi]
        x     = zeros(Float64, size(M,1))
        b     = zero(x)
        b[1:length(Fx)]                       .= Fx[:]
        b[length(Fx)+1:length(Fx)+length(Fy)] .= Fy[:]
        b[length(Fx)+length(Fy)+1:end]        .= Fp[:]
        epsi   = 1e-8
        noisy  = 2
        KSP_GCR_Stokes!( x, M, b, epsi, noisy, gamma, etac, Kuu, Kup, Kpu  )
        Vx[:,2:end-1] .+= x[NumVx]
        Vy[2:end-1,:] .+= x[NumVy]
        Pc            .+= x[NumP.+maximum(NumVy)]
    end
    return
end
export StokesSolver!

########

@views function DecoupledSolver!( Vx::Matrix{Float64}, Vy::Matrix{Float64}, Pc::Matrix{Float64}, NumVx::Matrix{Int64}, NumVy::Matrix{Int64}, NumP::Matrix{Int64}, fu::Vector{Float64}, fp::Vector{Float64}, Kuu::SparseMatrixCSC{Float64, Int64}, Kup::SparseMatrixCSC{Float64, Int64}, Kpu::SparseMatrixCSC{Float64, Int64}, etac::Matrix{Float64}, gamma::Float64 )
    # Decoupled solve
    ndofu = size(Kup,1)
    ndofp = size(Kup,2)
    @tturbo coef  = gamma*ones(length(etac))#.*etac[:]
    @tturbo Kppi  = spdiagm(coef)
    @tturbo Kuusc = Kuu - Kup*(Kppi*Kpu) # OK
    @tturbo PC    =  0.5*(Kuusc + Kuusc') 
    t = @elapsed Kf    = cholesky(Hermitian(PC), check = false)
    @printf("Cholesky took = %02.2e s\n", t)
    u     = zeros(Float64,ndofu)
    ru    = zeros(Float64,ndofu)
    fusc  = zeros(Float64,ndofu)
    p     = zeros(Float64,ndofp)
    rp    = zeros(Float64,ndofp)
    # Iterations
    for rit=1:10
        @tturbo ru   .= fu - Kuu*u - Kup*p;
        @tturbo rp   .= fp - Kpu*u;
        @printf("  --> Powell-Hestenes Iteration %02d\n  Momentum res.   = %2.2e\n  Continuity res. = %2.2e\n", rit, norm(ru)/sqrt(length(ru)), norm(rp)/sqrt(length(rp)))
        if norm(ru)/(length(ru)) < 1e-10 && norm(rp)/(length(ru)) < 1e-10
            break
        end
        @tturbo fusc .=  fu .- Kup*(Kppi*fp .+ p)
        @tturbo u    .= Kf\fusc
        @tturbo p   .+= Kppi*(fp .- Kpu*u)
    end
    Vx[:,2:end-1] .+= u[NumVx]
    Vy[2:end-1,:] .+= u[NumVy]
    Pc            .+= p[NumP]
end
export DecoupledSolver!

##############
function KSP_GCR_Stokes!( x::Vector{Float64}, M::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, eps::Float64, noisy::Int64, gamma::Float64, etac::Matrix{Float64}, Kuu::SparseMatrixCSC{Float64, Int64}, Kup::SparseMatrixCSC{Float64, Int64}, Kpu::SparseMatrixCSC{Float64, Int64} )
    # KSP GCR solver
    norm_r, norm0 = 0.0, 0.0
    N         = length(x)
    restart   = 25
    maxit     = 1000
    ncyc, its = 0, 0
    i1, i2, success=0,0,0
    # Arrays for coupled problem
    f      = zeros(Float64, N)
    v      = zeros(Float64, N)
    s      = zeros(Float64, N)
    val    = zeros(Float64, restart)
    VV     = zeros(Float64, (restart,N))
    SS     = zeros(Float64, (restart,N))
    # Coupled
    # Initial residual
    # Jx     = zeros(Float64, N)
    # Mx     = JacobianAction!(Jx, M, x; r,kv,T,fc,TW,TE,dx,n)
    f      = b - M*x 
    norm_r = norm(f)
    norm0  = norm_r;
    #
    ndofu = size(Kup,1)
    ndofp = size(Kup,2)
    coef  = gamma*ones(length(etac))#.*etac[:]
    Kppi  = spdiagm(coef)
    Kuusc = Kuu - Kup*(Kppi*Kpu) # OK
    PC    =  0.5*(Kuusc + Kuusc') 
    t = @elapsed Kf    = cholesky(Hermitian(PC),check = false)
    @printf("Cholesky took = %02.2e s\n", t)
    # Arrays for decoupled problem
    su    = zeros(Float64, ndofu)
    fusc  = zeros(Float64, ndofu)
    sp    = zeros(Float64, ndofp)
    fu    = zeros(Float64, ndofu)
    fp    = zeros(Float64, ndofp)
    @tturbo fu     .= f[1:ndofu]
    @tturbo fp     .= f[ndofu+1:end]
    if (noisy > 1) @printf("       %1.4d KSP GCR Residual %1.12e %1.12e\n", 0, norm_r, norm_r/norm0); end
    # Solving procedure
     while ( success == 0 && its<maxit ) 
        for i1=1:restart
            # Apply preconditioner, s = PC^{-1} f
            # s = PC\f
            @tturbo fusc .= fu  - Kup*(Kppi*fp + sp)
            @tturbo su   .= Kf\fusc
            @tturbo sp   .+= Kppi*(fp - Kpu*su)
            @tturbo s[1:ndofu]     .= su
            @tturbo s[ndofu+1:end] .= sp
            # Action of Jacobian on s: v = J*s
            # JacobianAction!(v, M, s; r,kv,T,fc,TW,TE,dx,n)
            @tturbo v .= M*s
            # Approximation of the Jv product
            for i2=1:i1
                @tturbo val[i2] = v' * VV[i2,:]
            end
            # Scaling
            for i2=1:i1
                @tturbo v .-= val[i2] * VV[i2,:]
                @tturbo s .-= val[i2] * SS[i2,:]
            end
            # -----------------
            @tturbo r_dot_v = f'*v
            nrm     = norm(v)
            r_dot_v = r_dot_v / nrm
            # -----------------
            fact    = 1.0/nrm
            @tturbo v     .*= fact
            @tturbo s     .*= fact
            # -----------------
            fact    = r_dot_v;
            @tturbo x     .+= fact*s
            @tturbo f     .-= fact*v
            # -----------------
            norm_r  = norm(f) 
            @tturbo fu     .= f[1:ndofu]
            @tturbo fp     .= f[ndofu+1:end]
            @printf("  --> Powell-Hestenes Iteration %02d\n  Momentum res.   = %2.2e\n  Continuity res. = %2.2e\n", its, norm(fu)/sqrt(length(fu)), norm(fp)/sqrt(length(fp)))
            if norm(fu)/(length(fu)) < 1e-10 && norm(fp)/(length(fu)) < 1e-10 #(norm_r < eps * norm0 )
                success = 1
                println("converged")
                break
            end
            # Store 
            @tturbo VV[i1,:] .= v
            @tturbo SS[i1,:] .= s
            its              += 1
        end
        its  += 1
        ncyc += 1
    end
    if (noisy>1) @printf("[%1.4d] %1.4d KSP GCR Residual %1.12e %1.12e\n", ncyc, its, norm_r, norm_r/norm0); end
    return its
end
export KSP_GCR_Stokes!

########