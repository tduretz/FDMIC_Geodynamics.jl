@views function StokesSolver!( f::Fields2D, params::ModelParameters, dom::ModelDomain, mat::MaterialParameters, BC::BoundaryConditions, Kuu::SparseMatrixCSC{Float64, Int64}, Kup::SparseMatrixCSC{Float64, Int64}, Kpu::SparseMatrixCSC{Float64, Int64}) #; dom=0, BC=0, phase_perc=0, nphase=0
    """ Solver """
    if params.JFNK==1 
        println("Automatically switching to KSP GCR solver since JFNK is requested")
        solver = 2 
    end
    # Get input residual
    fu = zeros(Float64, length(f.Fx) + length(f.Fy))
    fp = zeros(Float64, length(f.Fp))
    fu[1:length(f.Fx)]     .= f.Fx[:]
    fu[length(f.Fx)+1:end] .= f.Fy[:]
    fp .= f.Fp[:]
    #------------------------ PINNED PRESSURE
    if params.solver == -1
        # Pinned pressure case
        # Slightly compressible pressure block
        print("Pinned pressure")
        cP   = zeros(Float64,size(f.etac,1), size(f.etac,2))
        cP[1,1] = 1.0
        Kpu[1,:] .= 0.0
        f.Pc[1,1] = 0.0
        I    = f.NumP[:]
        J    = f.NumP[:]
        V    = cP[:]
        PP   = sparse( I, J, V)
        Kmat = [Kuu Kup; Kpu PP]
        F    = [fu; fp]
        dX   = Kmat\F
        println(size(f.Vx))
        println(size(f.NumVx))
        dVx = dX[ f.NumVx ]
        println( size(f.Vx[:,2:end-1]) )
        f.Vx[:,2:end-1] .= f.Vx[:,2:end-1] .+ dX[f.NumVx]
        f.Vy[2:end-1,:] .= f.Vy[2:end-1,:] .+ dX[f.NumVy] 
        f.Pc            .= f.Pc            .+ dX[f.NumP.+maximum(f.NumVy)]
        f.Pc            .-= mean(f.Pc)
    elseif params.solver == 0
        #------------------------ PINNED PRESSURE
        # Slightly compressible pressure block
        if params.comp==1
            print("The solver accounts for physical compressibility")
            cP   = 1.0./(f.Kc.*params.dt).*ones(Float64,size(f.etac,1), size(f.etac,2))
        else
            print("The solver does not account for physical compressibility (artificial only)")
            cP   = 1.0./params.gamma.*ones(Float64,size(f.etac,1), size(f.etac,2))
        end
        I    = f.NumP[:]
        J    = f.NumP[:]
        V    = cP[:]
        PP   = sparse( I, J, V)
        Kmat = [Kuu Kup; Kpu PP]
        F    = [fu; fp]
        dX   = Kmat\F
        println(size(f.Vx))
        println(size(f.NumVx))
        dVx = dX[ f.NumVx ]
        println( size(f.Vx[:,2:end-1]) )
        f.Vx[:,2:end-1] .= f.Vx[:,2:end-1] .+ dX[f.NumVx]
        f.Vy[2:end-1,:] .= f.Vy[2:end-1,:] .+ dX[f.NumVy] 
        f.Pc            .= f.Pc            .+ dX[f.NumP.+maximum(f.NumVy)]
        if params.comp==0 f.Pc            .-= mean(f.Pc) end
    elseif params.solver == 1
        DecoupledSolver!( f, params, fu, fp, Kuu, Kup, Kpu )
    elseif params.solver == 2 
        nVx = length(f.Fx)
        nVy = length(f.Fy)
        nP  = length(f.Fp)
        if params.comp==1
            @tturbo coef      = 1.0./(f.Kc[:].*params.dt).*ones(length(f.etac))#.*etac[:]
            @tturbo coef_inv  = f.Kc[:].*params.dt.*ones(length(f.etac))#.*etac[:]
        else
            @tturbo coef      =    0.0*ones(length(f.etac))#.*etac[:]
            @tturbo coef_inv  = params.gamma*ones(length(f.etac))#.*etac[:]
        end
        @tturbo Kpp   = spdiagm(coef)
        @tturbo Kppi  = spdiagm(coef_inv)
        @tturbo M     = [Kuu Kup; Kpu Kpp]
        dx    = zeros(Float64, size(M,1))
        x     = zero(dx)
        Vxi   = f.Vx[:,2:end-1]
        Vyi   = f.Vy[2:end-1,:]
        x[1:nVx]         .= Vxi[:]
        x[nVx+1:nVx+nVy] .= Vyi[:]
        x[nVx+nVy+1:end] .= f.Pc[:]
        b     = zero(dx)
        b[1:nVx]         .= f.Fx[:]
        b[nVx+1:nVx+nVy] .= f.Fy[:]
        b[nVx+nVy+1:end] .= f.Fp[:]
        epsi   = 1e-8
        noisy  = 2
        # ACHTUNG: here copy solution arrays since they are overwitten in the JFNK solver
        Vx, Vy, Pc = copy(f.Vx), copy(f.Vy), copy(f.Pc)
        if params.JFNK == 0
            KSP_GCR_Stokes!( dx, M, b, epsi, noisy, Kuu, Kup, Kpu, Kpp, Kppi )
        else
            KSP_GCR_Stokes2!( dx, M, x, b, epsi, noisy, Kuu, Kup, Kpu, Kpp, Kppi, params, dom, mat, f, BC )
        end
        # ACHTUNG: increment solutions from copied array
        f.Vx[:,2:end-1] .= Vx[:,2:end-1] .+ dx[f.NumVx]
        f.Vy[2:end-1,:] .= Vy[2:end-1,:] .+ dx[f.NumVy]
        f.Pc            .= Pc            .+ dx[f.NumP.+maximum(f.NumVy)]
    end
    return
end
export StokesSolver!

##############

@views function DecoupledSolver!( f::Fields2D, params::ModelParameters, fu::Vector{Float64}, fp::Vector{Float64}, Kuu::SparseMatrixCSC{Float64, Int64}, Kup::SparseMatrixCSC{Float64, Int64}, Kpu::SparseMatrixCSC{Float64, Int64} )
    # Decoupled solve
    ndofu = size(Kup,1)
    ndofp = size(Kup,2)
    if params.comp==1
        @tturbo coef      = 1.0./(f.Kc[:].*params.dt).*ones(length(f.etac))#.*etac[:]
        @tturbo coef_inv  = f.Kc[:].*params.dt.*ones(length(f.etac))#.*etac[:]
    else
        @tturbo coef      =    0.0*ones(length(f.etac))#.*etac[:]
        @tturbo coef_inv  = params.gamma*ones(length(f.etac))#.*etac[:]
    end
    @tturbo Kpp   = spdiagm(coef)
    @tturbo Kppi  = spdiagm(coef_inv)
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
        @tturbo ru   .= fu .- Kuu*u .- Kup*p;
        @tturbo rp   .= fp .- Kpu*u .- Kpp*p;
        @printf("  --> Powell-Hestenes Iteration %02d\n  Momentum res.   = %2.2e\n  Continuity res. = %2.2e\n", rit, norm(ru)/sqrt(length(ru)), norm(rp)/sqrt(length(rp)))
        if norm(ru)/(length(ru)) < 1e-13 && norm(rp)/(length(ru)) < 1e-13
            break
        end
        @tturbo fusc .=  fu .- Kup*(Kppi*fp .+ p)
        @tturbo u    .= Kf\fusc
        @tturbo p   .+= Kppi*(fp .- Kpu*u .- Kpp*p)
    end
    f.Vx[:,2:end-1] .= f.Vx[:,2:end-1] .+ u[f.NumVx]
    f.Vy[2:end-1,:] .= f.Vy[2:end-1,:] .+ u[f.NumVy]
    f.Pc            .= f.Pc            .+ p[f.NumP]
end
export DecoupledSolver!

##############

@views function KSP_GCR_Stokes!( x::Vector{Float64}, M::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, eps::Float64, noisy::Int64, Kuu::SparseMatrixCSC{Float64, Int64}, Kup::SparseMatrixCSC{Float64, Int64}, Kpu::SparseMatrixCSC{Float64, Int64}, Kpp::SparseMatrixCSC{Float64, Int64}, Kppi::SparseMatrixCSC{Float64, Int64} )
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
    f      = b - M*x 
    norm_r = norm(f)
    norm0  = norm_r;
    #
    ndofu = size(Kup,1)
    ndofp = size(Kup,2)
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

##############

function KSP_GCR_Stokes2!( x::Vector{Float64}, M::SparseMatrixCSC{Float64, Int64}, X::Vector{Float64}, b::Vector{Float64}, eps::Float64, noisy::Int64, Kuu::SparseMatrixCSC{Float64, Int64}, Kup::SparseMatrixCSC{Float64, Int64}, Kpu::SparseMatrixCSC{Float64, Int64}, Kpp::SparseMatrixCSC{Float64, Int64}, Kppi::SparseMatrixCSC{Float64, Int64}, params::ModelParameters, dom::ModelDomain, materials::MaterialParameters, f::Fields2D, BC::BoundaryConditions )
    # KSP GCR solver
    norm_r, norm0 = 0.0, 0.0
    N         = length(x)
    restart   = 30
    maxit     = 10#1000
    ncyc, its = 0, 0
    i1, i2, success=0,0,0
    epsi      = 1e-5
    ncx, ncy = dom.ncx, dom.ncy
    # Allocate tables
    Tau       = Tensor2D( zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+1, ncy+1), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0) ) 
    Eps       = Tensor2D( zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+1, ncy+1), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0) ) 
    F1        = zero(x)
    F2        = zero(x)
    # Arrays for coupled problem
    r      = zeros(Float64, N)
    v      = zeros(Float64, N)
    s      = zeros(Float64, N)
    val    = zeros(Float64, restart)
    VV     = zeros(Float64, (restart,N))
    SS     = zeros(Float64, (restart,N))
    # Coupled
    # Initial residual
    r      .= b# - Mx
    norm_r = norm(r)
    norm0  = norm_r;
    #
    ndofu = size(Kup,1)
    ndofp = size(Kup,2)
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
    @tturbo fu     .= r[1:ndofu]
    @tturbo fp     .= r[ndofu+1:end]
    if (noisy > 1) @printf("       %1.4d KSP GCR Residual %1.12e %1.12e\n", 0, norm_r, norm_r/norm0); end
    
    # res!( F1, X, f, Eps, Tau, dom, BC, params, materials, 1 ) 
    # println("|Fx| = ", norm(f.Fx)/length(f.Fx))
    # println("|Fy| = ", norm(f.Fy)/length(f.Fy))
    # println("|Fp| = ", norm(f.Fp)/length(f.Fp))

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
            # @tturbo v .= M*s

            # res!( F1, s, f, Eps, Tau, dom, BC, params, materials, 0 ) 
            # v .= F1

            res!( F1, X .+ epsi.*s, f, Eps, Tau, dom, BC, params, materials, 1 )
            res!( F2, X .- epsi.*s, f, Eps, Tau, dom, BC, params, materials, 1 )
            v .= (F1 .- F2) ./ ( 2epsi )
            
            # Approximation of the Jv product
            for i2=1:i1
                @tturbo val[i2] = v' * VV[i2,:]
            end
            # Scaling
            for i2=1:i1
                @tturbo v .-= val[i2] .* VV[i2,:]
                @tturbo s .-= val[i2] .* SS[i2,:]
            end
            # -----------------
            @tturbo r_dot_v = r'*v
            nrm      = sqrt(v'*v)
            r_dot_v /= nrm
            # -----------------
            @tturbo v     ./= nrm
            @tturbo s     ./= nrm
            # -----------------
            @tturbo x     .+= r_dot_v*s
            @tturbo r     .-= r_dot_v*v
            # -----------------
            norm_r  = sqrt(r'*r) 
            @tturbo fu     .= r[1:ndofu]
            @tturbo fp     .= r[ndofu+1:end]
            # @printf("  --> Powell-Hestenes Iteration %02d\n  Momentum res.   = %2.2e\n  Continuity res. = %2.2e\n", its, norm(fu)/sqrt(length(fu)), norm(fp)/sqrt(length(fp)))
            if norm(fu)/(length(fu)) < 1e-10 && norm(fp)/(length(fu)) < 1e-10 #(norm_r < eps * norm0 )
                success = 1
                # println("converged")
                break
            end
            # Store 
            @tturbo VV[i1,:] .= v[:]
            @tturbo SS[i1,:] .= s[:]
            its              += 1
        end
        its  += 1
        ncyc += 1
    end
    if (noisy>1) @printf("[%1.4d] %1.4d KSP GCR Residual %1.12e %1.12e\n", ncyc, its, norm_r, norm_r/norm0); end
    return its
end
export KSP_GCR_Stokes2!

##############

function Rheology_local!( f::Fields2D, Eps::Tensor2D, Tau::Tensor2D, params::ModelParameters, materials::MaterialParameters, BC::BoundaryConditions, ncx::Int64, ncy::Int64 )
    # Compute effective viscosity for each phase
    f.etac .= 0.0
    # f.Kc   .= 0.0
    for m=1:materials.nphase
        @tturbo eta    = 2.0* materials.eta0[m] * Eps.II.^(1.0/materials.n[m] - 1.0)
        @tturbo f.etac .+= f.phase_perc[m,:,:] .* eta
        # @tturbo f.Kc   .+= f.phase_perc[m,:,:] .* materials.K[m]
    end
    CentroidsToVertices!( f.etav, f.etac, ncx, ncy, BC )
    # println("min Eii: ", minimum(Eps.II), " --- max Eii: ", maximum(Eps.II))
    # println("min eta: ", minimum(f.etac), " --- max eta: ", maximum(f.etac))
end
export Rheology_local!

##############

function res!( r::Vector{Float64}, x::Vector{Float64}, f::Fields2D, Eps::Tensor2D, Tau::Tensor2D, dom::ModelDomain, BC::BoundaryConditions, params::ModelParameters, materials::MaterialParameters, nl::Int64 ) 
    f.Vx[:,2:end-1] .= x[f.NumVx]
    f.Vy[2:end-1,:] .= x[f.NumVy]
    f.Pc            .= x[f.NumP .+ maximum(f.NumVy)]
    SetBCs( f.Vx, f.Vy, BC )
    StrainRate!(Eps, f.Vx, f.Vy, dom)
    if nl==1
         Rheology_local!( f, Eps, Tau, params, materials, BC, dom.ncx, dom.ncy  )
    end
    Stress!( Tau, Eps, f )
    ResidualsComp!( f, Tau, Eps, BC, dom, params )
    nVx = length(f.Fx)
    nVy = length(f.Fy)
    r[1:nVx]         .= -f.Fx[:]
    r[nVx+1:nVx+nVy] .= -f.Fy[:]
    r[nVx+nVy+1:end] .= -f.Fp[:]
end
export res!