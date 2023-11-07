@views function StrainRate!( Eps::Tensor2D, Vx::Matrix{Float64}, Vy::Matrix{Float64}, dom::ModelDomain )
    dx, dy = dom.dx, dom.dy
     Eps.div    .= diff(Vx[:,2:end-1],dims=1)./dx .+ diff(Vy[2:end-1,:],dims=2)./dy 
    # ncx, ncy = size(div,1), size(div,2)
    #  for i=1:ncx
    #     for j=1:ncy
    #         div[i,j] = 1.0/dx * (Vx[i+1,j+1] - Vx[i,j+1]) + 1.0/dy * (Vy[i+1,j+1] - Vy[i+1,j])
    #     end
    # end
     Eps.xx .= diff(Vx[:,2:end-1],dims=1)./dx .- 1.0/3.0 .* Eps.div
     Eps.yy .= diff(Vy[2:end-1,:],dims=2)./dy .- 1.0/3.0 .* Eps.div
     Eps.zz .= .-(Eps.xx .+ Eps.yy)
     Eps.xy .= 0.5.*( diff(Vx,dims=2)./dy .+ diff(Vy,dims=1)/dx ) 
    VerticesToCentroids!( Eps.xy_c, Eps.xy )
     Eps.II .= sqrt.(0.5*(Eps.xx.^2 .+ Eps.xx.^2 .+ Eps.zz.^2) .+ Eps.xy_c.^2)
end
export StrainRate!

########

@views function Stress!( Tau::Tensor2D, Eps::Tensor2D, f::Fields2D )
     Tau.xx .= 2.0.*f.etac.*Eps.xx 
     Tau.yy .= 2.0.*f.etac.*Eps.yy
     Tau.zz .= 2.0.*f.etac.*Eps.zz  
     Tau.xy .= 2.0.*f.etav.*Eps.xy 
    VerticesToCentroids!( Tau.xy_c, Tau.xy )
     Tau.II .= sqrt.(0.5*(Tau.xx.^2 .+ Tau.xx.^2 .+ Tau.zz.^2) .+ Tau.xy_c.^2)
end
export Stress!

########

@views function StressVE!( Tau::Tensor2D, Tau0::Tensor2D, Eps::Tensor2D, f::Fields2D, dt::Float64 )
     Tau.xx .= 2f.etac.* ( Eps.xx .+ Tau0.xx ./ ( 2f.Gc.*dt ) ) 
     Tau.yy .= 2f.etac.* ( Eps.yy .+ Tau0.yy ./ ( 2f.Gc.*dt ) ) 
     Tau.zz .= 2f.etac.* ( Eps.zz .+ Tau0.zz ./ ( 2f.Gc.*dt ) ) 
     Tau.xy .= 2f.etav.* ( Eps.xy .+ Tau0.xy ./ ( 2f.Gv.*dt ) ) 
    VerticesToCentroids!( Tau.xy_c, Tau.xy )
     Tau.II .= sqrt.(0.5*((f.Damc.*Tau.xx).^2 .+ (f.Damc.*Tau.yy).^2 .+ (f.Damc.*Tau.zz).^2) .+ (f.Damc.*Tau.xy_c).^2)
end
export StressVE!

########
function StrainEnergy!( We::Matrix{Float64}, Tau::Tensor2D, Eps::Tensor2D, Str0::Tensor2D, Str::Tensor2D, P::Matrix{Float64}, dt::Float64, f::Fields2D, dx::Float64, dy::Float64 )
    Str.xx   .= Str0.xx   .+ dt.*(Eps.xx .+ 1.0/3.0.*Eps.div)
    Str.yy   .= Str0.yy   .+ dt.*(Eps.yy .+ 1.0/3.0.*Eps.div)
    Str.zz   .= Str0.zz   .+ dt.*(Eps.zz .+ 1.0/3.0.*Eps.div)
    Str.xy_c .= Str0.xy_c .+ dt.*(Eps.xy_c)
     Str.II .= sqrt.(0.5*(Str.xx.^2 .+ Str.xx.^2 .+ Str.zz.^2) .+ Str.xy_c.^2)
    We  .= (Tau.xx .- P) .* Str.xx .+ (Tau.yy .- P) .* Str.yy .+ (Tau.zz .- P) .* Str.zz .+ 2.0.*Tau.xy_c .* Str.xy_c
    We .*= 0.5  
    # We ./= f.Damc
    We[f.We0.>We] .= f.We0[f.We0.>We]

    # Exx = diff(f.Vx[:,2:end-1],dims=1)./dx * dt
    # Eyy = diff(f.Vy[2:end-1,:],dims=2)./dy * dt
    # Exyv = 0.5.*( diff(f.Vx,dims=2)./dy .+ diff(f.Vy,dims=1)/dx ) * dt  
    # Exyc = zeros(Float64, size(Exx))
    # VerticesToCentroids!( Exyc, Exyv )

    # Sxx = (f.Kc.+4.0./3.0*f.Gc).*Exx + (f.Kc.-2.0./3.0*f.Gc).*Eyy
    # Syy = (f.Kc.-2.0./3.0*f.Gc).*Exx + (f.Kc.+4.0./3.0*f.Gc).*Eyy
    # Szz = (f.Kc.-4.0./3.0*f.Gc).*Exx + (f.Kc.-2.0./3.0*f.Gc).*Eyy
    # Sxy = 2.0*f.Gc.*Exyc

    # We .= 2.0*Sxy.*Exyc .+ Sxx.*Exx .+ Syy.*Eyy
end
export StrainEnergy!
########

@views function Residuals!( f::Fields2D, Tau::Tensor2D, Eps::Tensor2D, BC::BoundaryConditions, dom::ModelDomain, params::ModelParameters )
    dx, dy = dom.dx, dom.dy
     f.Fx .= 0.0
     f.Fy .= 0.0
     f.Fp .= 0.0
    if BC.periodix==0
        f.Fx[2:end-1,:] .= ( diff(Tau.xx .- f.Pc, dims=1)./dx .+ diff(Tau.xy[2:end-1,:], dims=2)/dy )
    else
        Sxx_ex             = zeros(ncx+1,ncy)
         Sxx_ex[2:end-0,:] .= -f.Pc .+ Tau.xx
        Sxx_ex[      1,:] .= -f.Pc[end,:] .+ Tau.xx[end,:]
        #Sxx_ex[    end,:] .= -Pc[  1,:] .+ Tau.xx[  1,:] # Do not assemble last column
        f.Fx .= ( diff(Sxx_ex, dims=1)./dx .+ diff(Tau.xy[1:end-1,:], dims=2)/dy )
    end
    f.Fy[:,2:end-1] .= ( diff(Tau.yy .- f.Pc, dims=2)./dy .+ diff(Tau.xy[:,2:end-1], dims=1)/dx )
     f.Fp            .= -Eps.div
    # # For periodic
    # if BC.periodix==1
    #     Fx = Fx[1:end-1,:]
    # end
end
export Residuals!

########

@views function ResidualsComp!( f::Fields2D, Tau::Tensor2D, Eps::Tensor2D, BC::BoundaryConditions, dom::ModelDomain, params::ModelParameters )
    dx, dy, ncx, ncy = dom.dx, dom.dy, dom.ncx, dom.ncy 
     f.Fx .= 0.0
     f.Fy .= 0.0
     f.Fp .= 0.0
    if BC.periodix==0
        f.Fx[2:end-1,:] .= ( diff(Tau.xx .- f.Pc, dims=1)./dx .+ diff(Tau.xy[2:end-1,:], dims=2)/dy )
    else
        Sxx_ex             = zeros(ncx+1,ncy)
         Sxx_ex[2:end-0,:] .= -f.Pc .+ Tau.xx
        Sxx_ex[      1,:] .= -f.Pc[end,:] .+ Tau.xx[end,:]
        #Sxx_ex[    end,:] .= -Pc[  1,:] .+ Tau.xx[  1,:] # Do not assemble last column
        f.Fx .= ( diff(Sxx_ex, dims=1)./dx .+ diff(Tau.xy[1:end-1,:], dims=2)/dy )
    end
    f.Fy[:,2:end-1] .= ( diff(Tau.yy .- f.Pc, dims=2)./dy .+ diff(Tau.xy[:,2:end-1], dims=1)/dx ) .- params.gy.*0.5.*(f.rhoc[:,1:end-1] .+ f.rhoc[:,2:end])
    if params.comp == 1
         f.Fp .= -Eps.div .- (f.Pc .- f.Pc0) ./ (f.Kc.*params.dt)
    else
         f.Fp .= -Eps.div
    end
end
export ResidualsComp!

########

@views function ResidualsCompDam!( f::Fields2D, Tau::Tensor2D, Eps::Tensor2D, BC::BoundaryConditions, dom::ModelDomain, params::ModelParameters )
    dx, dy, ncx, ncy = dom.dx, dom.dy, dom.ncx, dom.ncy 
     f.Fx .= 0.0
     f.Fy .= 0.0
     f.Fp .= 0.0
    if BC.periodix==0
        f.Fx[2:end-1,:] .= ( diff(f.Damc.*(Tau.xx .- f.Pc), dims=1)./dx .+ diff(f.Damv[2:end-1,:].*Tau.xy[2:end-1,:], dims=2)/dy )
    else
        Sxx_ex             = zeros(ncx+1,ncy)
         Sxx_ex[2:end-0,:] .= -f.Pc .+ Tau.xx
        Sxx_ex[      1,:] .= -f.Pc[end,:] .+ Tau.xx[end,:]
        #Sxx_ex[    end,:] .= -Pc[  1,:] .+ Tau.xx[  1,:] # Do not assemble last column
        f.Fx .= ( diff(Sxx_ex, dims=1)./dx .+ diff(Tau.xy[1:end-1,:], dims=2)/dy )
    end
    rhoy = 0.5*( f.rhoc[:,1:end-1] .+ f.rhoc[:,2:end-0] )
    f.Fy[:,2:end-1] .=  params.gy.*rhoy .+ ( diff(f.Damc.*(Tau.yy .- f.Pc), dims=2)./dy .+ diff(f.Damv[:,2:end-1].*Tau.xy[:,2:end-1], dims=1)/dx )
    println(mean(rhoy))
    println(params.gy)
    if params.comp == 1
         f.Fp .= -Eps.div .- (f.Pc .- f.Pc0) ./ (f.Kc.*params.dt)
    else
         f.Fp .= -Eps.div
    end
end
export ResidualsCompDam!

########

@views function SetBCs( Vx::Matrix{Float64}, Vy::Matrix{Float64}, BC::BoundaryConditions )
    # Boundaries
    if BC.Vx.type_S == 22 Vx[:,  1] .= Vx[:,    2] end
    if BC.Vx.type_S == 11 Vx[:,  1] .= 2.0.*BC.Vx.Dir_S .- Vx[:,    2] end
    if BC.Vx.type_N == 22 Vx[:,end] .= Vx[:,end-1] end
    if BC.Vx.type_N == 11 Vx[:,end] .= 2.0.*BC.Vx.Dir_N .- Vx[:,end-1] end
    if BC.periodix == 1
        Vy[1,  :] .= Vy[end-1,:]
        Vy[end,:] .= Vy[    2,:]
    else
        if BC.Vy.type_W == 22  Vy[1,  :] .= Vy[2,    :] end
        if BC.Vy.type_W == 11  Vy[1,  :] .= 2.0.*BC.Vy.Dir_W .- Vy[2,    :] end
        if BC.Vy.type_E == 22  Vy[end,:] .= Vy[end-1,:] end
        if BC.Vy.type_E == 11  Vy[end,:] .= 2.0.*BC.Vy.Dir_E .- Vy[end-1,:] end
    end
end
export SetBCs

########

@views function NumberingStokes!( f::Fields2D, BC::BoundaryConditions, dom::ModelDomain )
    ncx, ncy = dom.ncx, dom.ncy
    # Numbering
    if BC.periodix==0
        f.NumVx     = collect(reshape(1:(ncx+1)*ncy,ncx+1,ncy))
    else
        f.NumVx             = zeros(Int64,(ncx+1),ncy)
        f.NumVx[1:end-1,:] .= collect(reshape(1:(ncx+0)*ncy,ncx,ncy))
        f.NumVx[end,:]     .= f.NumVx[1,:]
    end
    f.NumVy     = collect(reshape(1:ncx*(ncy+1),ncx,ncy+1) .+ maximum(f.NumVx))
    f.NumP      = collect(reshape(1:(ncx)*ncy,ncx,ncy))
    return
end
export NumberingStokes!

########

@views function StokesAssembly( BC::BoundaryConditions, f::Fields2D, DirScale::Float64, dom::ModelDomain )

    ncx, ncy = dom.ncx, dom.ncy
    dx, dy = dom.dx, dom.dy
    NumVx, NumVy, NumP = f.NumVx, f.NumVy, f.NumP
    etac, etav = f.Damc.*f.etac, f.Damv.*f.etav

    if BC.periodix==0
        nxvx = ncx+1
        sx   = 0
    else
        nxvx = ncx+1#ncx # remove last column
        sx   = 0#1   # shift 
    end

    # Connectivity
    iVxC      =  ones(Int64, nxvx, ncy);
    iVxC     .=  NumVx[1:end-sx,:]
    iVxW      =  ones(Int64, nxvx, ncy); iVxW[2:end-0,: ] = NumVx[1:end-1,:]
    iVxE      =  ones(Int64, nxvx, ncy); iVxE[1:end-1,: ] = NumVx[2:end-0,:]        
    iVxS      =  ones(Int64, nxvx, ncy); iVxS[: ,2:end-0] = NumVx[:,1:end-1]
    iVxN      =  ones(Int64, nxvx, ncy); iVxN[: ,1:end-1] = NumVx[:,2:end-0]
    iVySW     =  ones(Int64, nxvx, ncy); iVySW[2:end-0,:] = NumVy[:,1:end-1]
    iVySE     =  ones(Int64, nxvx, ncy); iVySE[1:end-1,:] = NumVy[:,1:end-1]
    iVyNW     =  ones(Int64, nxvx, ncy); iVyNW[2:end-0,:] = NumVy[:,2:end-0]
    iVyNE     =  ones(Int64, nxvx, ncy); iVyNE[1:end-1,:] = NumVy[:,2:end-0]
    iPW       =  ones(Int64, nxvx, ncy); iPW[2:end-0,:]   = NumP[:,:]
    iPE       =  ones(Int64, nxvx, ncy); iPE[1:end-1,:]   = NumP[:,:]

    # Viscosity coefficients
    etaW      = zeros(Float64, nxvx, ncy); etaW[2:end-0,:] = etac[1:end-0,:]
    etaE      = zeros(Float64, nxvx, ncy); etaE[1:end-1,:] = etac[1:end-0,:]
    etaS      = zeros(Float64, nxvx, ncy); etaS[:,1:end-0] = etav[:,1:end-1] 
    etaN      = zeros(Float64, nxvx, ncy); etaN[:,1:end-0] = etav[:,2:end-0]

    DamW =  ones(Float64, size(NumVx)); DamW[2:end-0,:]   .= f.Damc
    DamE =  ones(Float64, size(NumVx)); DamE[1:end-1,:]   .= f.Damc

    if BC.periodix==1
        DamW[  1,:] = f.Damc[end,:]
        DamE[end,:] = f.Damc[1,:]
        etaW[  1,:] = etac[end,:]
        etaE[end,:] = etac[1,:]
        iVxW[  1,:] = NumVx[end-1,:]
        iVxE[end,:] = NumVx[    2,:]
        iVySW[  1,:]  = NumVy[end,1:end-1]
        iVyNW[  1,:]  = NumVy[end,2:end-0]
        iVySE[end,:]  = NumVy[  1,1:end-1]
        iVyNE[end,:]  = NumVy[  1,2:end-0]
        iPW[  1,:] = NumP[end,:]
        iPE[end,:] = NumP[  1,:]
    end

    # Finite difference coefficients
     cVxC  = -(-1.0.*etaN./dy - 1.0.*etaS./dy)./dy - (-4/3*etaE./dx - 4/3*etaW./dx)./dx
     cVxW  = -4/3*etaW./dx.^2
     cVxE  = -4/3*etaE./dx.^2
     cVxS  = -1.0*etaS./dy.^2
     cVxN  = -1.0*etaN./dy.^2
     cVySW = -1.0*etaS./(dx.*dy) + 2/3*etaW./(dx.*dy)
     cVySE = -2/3*etaE./(dx.*dy) + 1.0*etaS./(dx.*dy)
     cVyNW = 1.0*etaN./(dx.*dy) - 2/3*etaW./(dx.*dy)
     cVyNE = 2/3*etaE./(dx.*dy) - 1.0*etaN./(dx.*dy)
     cPW   = -DamW./dx .*  ones(Float64, nxvx, ncy)
     cPE   =  DamE./dx .*  ones(Float64, nxvx, ncy)

    if BC.Vx.type_S==11
        cVxC[:,  1] .-= cVxS[:,  1]
    end
    if BC.Vx.type_N==11
        cVxC[:,end] .-= cVxN[:,end]
    end
    if BC.Vx.type_S==22
        cVxC[:,  1] .+= cVxS[:,  1]
    end
    if BC.Vx.type_N==22
        cVxC[:,end] .+= cVxN[:,end]
    end
    cVxS[:,  1] .= 0.0
    cVxN[:,end] .= 0.0

    # Symmetry - kill Dirichlet connections
    cVySW[  :,  1] .= 0.0
    cVySE[  :,  1] .= 0.0
    cVyNW[  :,end] .= 0.0
    cVyNE[  :,end] .= 0.0

    if BC.periodix==0
        cVxS[:,  1] .= 0.0; 
        cVxN[:,end] .= 0.0; 
        cVxC[1,:]  .= DirScale; cVxC[end,:]  .= DirScale
        cVxW[1,:]  .= 0.0;      cVxW[end,:]  .= 0.0
        cVxE[1,:]  .= 0.0;      cVxE[end,:]  .= 0.0
        cVxS[1,:]  .= 0.0;      cVxS[end,:]  .= 0.0
        cVxN[1,:]  .= 0.0;      cVxN[end,:]  .= 0.0
        cVySW[1,:] .= 0.0;      cVySW[end,:] .= 0.0
        cVyNW[1,:] .= 0.0;      cVyNW[end,:] .= 0.0
        cVySE[1,:] .= 0.0;      cVySE[end,:] .= 0.0
        cVyNE[1,:] .= 0.0;      cVyNE[end,:] .= 0.0
        cPW[1,:]   .= 0.0;      cPW[end,:]   .= 0.0
        cPE[1,:]   .= 0.0;      cPE[end,:]   .= 0.0
        # Symmetry - kill Dirichlet connections
        cVxW[    2,:] .= 0.0
        cVxE[end-1,:] .= 0.0
    end

    ###################

    # Connectivity
    iVyC      = NumVy
    iVyW      =  ones(Int64, size(NumVy)); iVyW[2:end-0,: ] = NumVy[1:end-1,:]
    iVyE      =  ones(Int64, size(NumVy)); iVyE[1:end-1,: ] = NumVy[2:end-0,:]
    iVyS      =  ones(Int64, size(NumVy)); iVyS[: ,2:end-0] = NumVy[:,1:end-1]
    iVyN      =  ones(Int64, size(NumVy)); iVyN[: ,1:end-1] = NumVy[:,2:end-0]
    iVxSW     =  ones(Int64, size(NumVy)); iVxSW[:,2:end-0] = NumVx[1:end-1,:]
    iVxSE     =  ones(Int64, size(NumVy)); iVxSE[:,2:end-0] = NumVx[2:end-0,:]
    iVxNW     =  ones(Int64, size(NumVy)); iVxNW[:,1:end-1] = NumVx[1:end-1,:]
    iVxNE     =  ones(Int64, size(NumVy)); iVxNE[:,1:end-1] = NumVx[2:end-0,:]
    iPS       =  ones(Int64, size(NumVy)); iPS[:,2:end-0]   = NumP
    iPN       =  ones(Int64, size(NumVy)); iPN[:,1:end-1]   = NumP

    if BC.periodix==1
        iVyW[  1,:] = NumVy[end,:]
        iVyE[end,:] = NumVy[  1,:]
    end

    # Viscosity coefficients
    etaS      = zeros(size(NumVy)); etaS[:,2:end-0] = etac[:,1:end-0]
    etaN      = zeros(size(NumVy)); etaN[:,1:end-1] = etac[:,1:end-0]
    etaW      = zeros(size(NumVy)); etaW[1:end-0,:] = etav[1:end-1,:] 
    etaE      = zeros(size(NumVy)); etaE[1:end-0,:] = etav[2:end-0,:]
    # Finite difference coefficients
     cVyC  = -(-4/3*etaN./dy - 4/3*etaS./dy)./dy - (-1.0*etaE./dx - 1.0*etaW./dx)./dx
     cVyW  = -1.0*etaW./dx.^2
     cVyE  = -1.0*etaE./dx.^2
     cVyS  = -4/3*etaS./dy.^2
     cVyN  = -4/3*etaN./dy.^2
     cVxSW = 2/3*etaS./(dx.*dy) - 1.0*etaW./(dx.*dy)
     cVxSE = 1.0*etaE./(dx.*dy) - 2/3*etaS./(dx.*dy)
     cVxNW = -2/3*etaN./(dx.*dy) + 1.0*etaW./(dx.*dy)
     cVxNE = -1.0*etaE./(dx.*dy) + 2/3*etaN./(dx.*dy)
    DamS =  ones(Float64, size(NumVy)); DamS[:,2:end-0]   .= f.Damc
    DamN =  ones(Float64, size(NumVy)); DamN[:,1:end-1]   .= f.Damc
     cPS   = -DamS./dy .* ones(size(NumVy)); cPS[:,  1] .= 0.0;  cPS[:,end] .= 0.0
     cPN   =  DamN./dy .* ones(size(NumVy)); cPN[:,  1] .= 0.0;  cPN[:,end] .= 0.0

    if BC.periodix==0
        if BC.Vy.type_W==11
            cVyC[  1,:] .-= cVyW[  1,:]
        end
        if BC.Vy.type_W==22
            cVyC[  1,:] .+= cVyW[  1,:]
        end
        if BC.Vy.type_E==11
            cVyC[end,:] .-= cVyE[end,:]
        end
        if BC.Vy.type_E==22
            cVyC[end,:] .+= cVyE[end,:]
        end
        cVyW[  1,:] .= 0.0
        cVyE[end,:] .= 0.0
    end

    # N-S Dirichlet nodes
    cVyC[:,1]  .= DirScale; cVyC[:,end]  .= DirScale
    cVyW[:,1]  .= 0.0;      cVyW[:,end]  .= 0.0
    cVyE[:,1]  .= 0.0;      cVyE[:,end]  .= 0.0
    cVyS[:,1]  .= 0.0;      cVyS[:,end]  .= 0.0
    cVyN[:,1]  .= 0.0;      cVyN[:,end]  .= 0.0
    cVxSW[:,1] .= 0.0;      cVxSW[:,end] .= 0.0
    cVxNW[:,1] .= 0.0;      cVxNW[:,end] .= 0.0
    cVxSE[:,1] .= 0.0;      cVxSE[:,end] .= 0.0
    cVxNE[:,1] .= 0.0;      cVxNE[:,end] .= 0.0
    # Symmetry - kill Dirichlet connections
    cVyS[:,     2] .= 0.0
    cVyN[:, end-1] .= 0.0

    if BC.periodix==0
        # cVyW[  1,:] .= 0.0;
        # cVyE[end,:] .= 0.0;

        # cVyC[:,1] .= 1e3; cVyC[:,end] .= 1e3
        # cVyW[:,1] .= 0.0; cVyW[:,end] .= 0.0
        # cVyE[:,1] .= 0.0; cVyE[:,end] .= 0.0
        # cVyS[:,1] .= 0.0; cVyS[:,end] .= 0.0
        # cVyN[:,1] .= 0.0; cVyN[:,end] .= 0.0
        # cVxSW[:,1] .= 0.0; cVxSW[:,end] .= 0.0
        # cVxNW[:,1] .= 0.0; cVxNW[:,end] .= 0.0
        # cVxSE[:,1] .= 0.0; cVxSE[:,end] .= 0.0
        # cVxNE[:,1] .= 0.0; cVxNE[:,end] .= 0.0
        # # Symmetry - kill Dirichlet connections
        # cVyS[:,   2]   .= 0.0
        # cVyN[:,end-1]  .= 0.0
        # cVxSW[  1,  :] .= 0.0
        # cVxSE[end,  :] .= 0.0
        # cVxNW[  1,  :] .= 0.0
        # cVxNE[end,  :] .= 0.0
        cVyW[  1,:] .= 0.0;
        cVyE[end,:] .= 0.0;
        # Symmetry - kill Dirichlet connections
        cVxSW[  1,  :] .= 0.0
        cVxSE[end,  :] .= 0.0
        cVxNW[  1,  :] .= 0.0
        cVxNE[end,  :] .= 0.0
    end

    if BC.periodix==1
        # Remove redundant Vx equation on the right side to make Kuu matrix symmetric positive definite
        iVxC  = collect(iVxC[1:end-1,:]);  cVxC  = collect(cVxC[1:end-1,:])
        iVxW  = collect(iVxW[1:end-1,:]);  cVxW  = collect(cVxW[1:end-1,:])
        iVxE  = collect(iVxE[1:end-1,:]);  cVxE  = collect(cVxE[1:end-1,:])
        iVxS  = collect(iVxS[1:end-1,:]);  cVxS  = collect(cVxS[1:end-1,:])
        iVxN  = collect(iVxN[1:end-1,:]);  cVxN  = collect(cVxN[1:end-1,:])
        iVySW = collect(iVySW[1:end-1,:]); cVySW = collect(cVySW[1:end-1,:])
        iVySE = collect(iVySE[1:end-1,:]); cVySE = collect(cVySE[1:end-1,:])
        iVyNW = collect(iVyNW[1:end-1,:]); cVyNW = collect(cVyNW[1:end-1,:])
        iVyNE = collect(iVyNE[1:end-1,:]); cVyNE = collect(cVyNE[1:end-1,:])
        iPW   = collect(iPW[1:end-1,:]);   cPW   = collect(cPW[1:end-1,:])
        iPE   = collect(iPE[1:end-1,:]);   cPE   = collect(cPE[1:end-1,:])
        # ncx = size(cVxC,1)-1
        # ncy = size(cVxC,2)
        # iVxC1  = zeros(Int64,(ncx,ncy)); cVxC1  = zeros(Float64,(ncx,ncy))
        # iVxW1  = zeros(Int64,(ncx,ncy)); cVxW1  = zeros(Float64,(ncx,ncy))
        # iVxE1  = zeros(Int64,(ncx,ncy)); cVxE1  = zeros(Float64,(ncx,ncy))
        # iVxS1  = zeros(Int64,(ncx,ncy)); cVxS1  = zeros(Float64,(ncx,ncy))
        # iVxN1  = zeros(Int64,(ncx,ncy)); cVxN1  = zeros(Float64,(ncx,ncy));
        # iVySW1  = zeros(Int64,(ncx,ncy)); cVySW1  = zeros(Float64,(ncx,ncy));
        # iVySE1  = zeros(Int64,(ncx,ncy)); cVySE1  = zeros(Float64,(ncx,ncy));
        # iVyNW1  = zeros(Int64,(ncx,ncy)); cVyNW1  = zeros(Float64,(ncx,ncy));
        # iVyNE1  = zeros(Int64,(ncx,ncy)); cVyNE1  = zeros(Float64,(ncx,ncy));
        # iPW1  = zeros(Int64,(ncx,ncy)); cPW1  = zeros(Float64,(ncx,ncy));
        # iPE1  = zeros(Int64,(ncx,ncy)); cPE1  = zeros(Float64,(ncx,ncy));
        # iVxC1  .= iVxC[1:end-1,:];  cVxC1  .= cVxC[1:end-1,:]
        # iVxW1  .= iVxW[1:end-1,:];  cVxW1  .= cVxW[1:end-1,:]
        # iVxE1  .= iVxE[1:end-1,:];  cVxE1  .= cVxE[1:end-1,:]
        # iVxS1  .= iVxS[1:end-1,:];  cVxS1  .= cVxS[1:end-1,:]
        # iVxN1  .= iVxN[1:end-1,:];  cVxN1  .= cVxN[1:end-1,:]
        # iVySW1 .= iVySW[1:end-1,:]; cVySW1 .= cVySW[1:end-1,:]
        # iVySE1 .= iVySE[1:end-1,:]; cVySE1 .= cVySE[1:end-1,:]
        # iVyNW1 .= iVyNW[1:end-1,:]; cVyNW1 .= cVyNW[1:end-1,:]
        # iVyNE1 .= iVyNE[1:end-1,:]; cVyNE1 .= cVyNE[1:end-1,:]
        # iPW1   .= iPW[1:end-1,:];   cPW1   .= cPW[1:end-1,:]
        # iPE1   .= iPE[1:end-1,:];   cPE1   .= cPE[1:end-1,:]
    end

    # Sparse matrix Kuu
    nVx = size(cVxC,1)*size(cVxC,2)
    nVy = size(cVyC,1)*size(cVyC,2)
    I   = zeros(  Int64, 9*(nVx + nVy) )
    J   = zeros(  Int64, 9*(nVx + nVy) )
    V   = zeros(Float64, 9*(nVx + nVy) )

    #------------------- Vx
    FillCoefficients!(I, J, V, 0*nVx, iVxC[:], iVxC[:], cVxC[:])
    FillCoefficients!(I, J, V, 1*nVx, iVxC[:], iVxW[:], cVxW[:])
    FillCoefficients!(I, J, V, 2*nVx, iVxC[:], iVxE[:], cVxE[:])
    FillCoefficients!(I, J, V, 3*nVx, iVxC[:], iVxS[:], cVxS[:])
    FillCoefficients!(I, J, V, 4*nVx, iVxC[:], iVxN[:], cVxN[:])
    FillCoefficients!(I, J, V, 5*nVx, iVxC[:], iVySW[:], cVySW[:])
    FillCoefficients!(I, J, V, 6*nVx, iVxC[:], iVySE[:], cVySE[:])
    FillCoefficients!(I, J, V, 7*nVx, iVxC[:], iVyNW[:], cVyNW[:])
    FillCoefficients!(I, J, V, 8*nVx, iVxC[:], iVyNE[:], cVyNE[:])
    #------------------- Vy
    FillCoefficients!(I, J, V, 9*nVx+0*nVy, iVyC[:], iVyC[:], cVyC[:])
    FillCoefficients!(I, J, V, 9*nVx+1*nVy, iVyC[:], iVyW[:], cVyW[:])
    FillCoefficients!(I, J, V, 9*nVx+2*nVy, iVyC[:], iVyE[:], cVyE[:])
    FillCoefficients!(I, J, V, 9*nVx+3*nVy, iVyC[:], iVyS[:], cVyS[:])
    FillCoefficients!(I, J, V, 9*nVx+4*nVy, iVyC[:], iVyN[:], cVyN[:])
    FillCoefficients!(I, J, V, 9*nVx+5*nVy, iVyC[:], iVxSW[:], cVxSW[:])
    FillCoefficients!(I, J, V, 9*nVx+6*nVy, iVyC[:], iVxSE[:], cVxSE[:])
    FillCoefficients!(I, J, V, 9*nVx+7*nVy, iVyC[:], iVxNW[:], cVxNW[:])
    FillCoefficients!(I, J, V, 9*nVx+8*nVy, iVyC[:], iVxNE[:], cVxNE[:])
    #------------------- Assemble
    Kuu = sparse( I, J, V)
    droptol!(Kuu, 1e-9)
    # display(UnicodePlots.spy(K))

    # Sparse matrix Kup
    I   = zeros(  Int64, 2*(nVx + nVy) )
    J   = zeros(  Int64, 2*(nVx + nVy) )
    V   = zeros(Float64, 2*(nVx + nVy) )
    FillCoefficients!(I, J, V, 0*nVx+0*nVy, iVxC[:], iPW[:], cPW[:])
    FillCoefficients!(I, J, V, 1*nVx+0*nVy, iVxC[:], iPE[:], cPE[:])
    FillCoefficients!(I, J, V, 2*nVx+0*nVy, iVyC[:], iPS[:], cPS[:])
    FillCoefficients!(I, J, V, 2*nVx+1*nVy, iVyC[:], iPN[:], cPN[:])
    Kup = sparse( I, J, V)
    droptol!(Kup, 1e-9)

    # Sparse matrix Kpu: Coefficients
    iVxW  =  ones(  Int64, size(NumP)); iVxW .= NumVx[1:end-1,:]
    iVxE  =  ones(  Int64, size(NumP)); iVxE .= NumVx[2:end-0,:]
    iVyS  =  ones(  Int64, size(NumP)); iVyS .= NumVy[:,1:end-1]
    iVyN  =  ones(  Int64, size(NumP)); iVyN .= NumVy[:,2:end-0]
    cVxW  = zeros(Float64, size(NumP)); cVxW .= -1.0 ./ dx
    cVxE  = zeros(Float64, size(NumP)); cVxE .=  1.0 ./ dx
    cVyS  = zeros(Float64, size(NumP)); cVyS .= -1.0 ./ dy
    cVyN  = zeros(Float64, size(NumP)); cVyN .=  1.0 ./ dy

    # Kill Dirichlet connections
    cVyS[:,1] .= 0.0; cVyN[:,end] .= 0.0
    if BC.periodix==0
        cVxW[1,:] .= 0.0; cVxE[end,:] .= 0.0
    end
    # Sparse matrix Kpu: Triplets
    nP = length(NumP)
    I  = zeros(  Int64, 4*nP)
    J  = zeros(  Int64, 4*nP) 
    V  = zeros(Float64, 4*nP) 
    FillCoefficients!(I, J, V, 0*nP, NumP[:], iVxW[:], cVxW[:])
    FillCoefficients!(I, J, V, 1*nP, NumP[:], iVxE[:], cVxE[:])
    FillCoefficients!(I, J, V, 2*nP, NumP[:], iVyS[:], cVyS[:])
    FillCoefficients!(I, J, V, 3*nP, NumP[:], iVyN[:], cVyN[:])
    Kpu = sparse( I, J, V)
    droptol!(Kpu, 1e-9)
    # display(UnicodePlots.spy(Div))
    return Kuu, Kup, Kpu 
end
export StokesAssembly

########

@views function SetInitialVelocity!( Vx::Matrix{Float64}, Vy::Matrix{Float64}, BC::BoundaryConditions, xv::LinRange{Float64}, yv::LinRange{Float64}, xce::LinRange{Float64}, yce::LinRange{Float64}, xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64, ncx::Int64, ncy::Int64 )

    # Initial guess - Vx
    Vx[1    ,:] .= BC.Vx.Dir_W
    Vx[end  ,:] .= BC.Vx.Dir_E
    Vx[:,    2] .= BC.Vx.Dir_S
    Vx[:,end-1] .= BC.Vx.Dir_N
    Vy[2    ,:] .= BC.Vy.Dir_W
    Vy[end-1,:] .= BC.Vy.Dir_E
    Vy[:,    1] .= BC.Vy.Dir_S
    Vy[:,  end] .= BC.Vy.Dir_N

    Vx_ini = zeros(ncx+1,ncy+2)
    Vx_ini .= Vx
    Vy_ini = zeros(ncx+2,ncy+1)
    Vy_ini .= Vy

    for i=1:size(Vx,1)
        for j=1:size(Vx,2)
            wW = 1.0 - (xv[i]-xmin)/(xmax-xmin)
            Vx[i,j] = wW * BC.Vx.Dir_W[j] + (1.0-wW) * BC.Vx.Dir_E[j]
            if i>1 && i<size(Vx,1)
                wS = 1.0 - (yce[j]-ymin)/(ymax-ymin)
                Vx[i,j] += wS * BC.Vx.Dir_S[i] + (1.0-wS) * BC.Vx.Dir_N[i]
            end
        end
    end

    # Initial guess - Vy
    for i=1:size(Vy,1)
        for j=1:size(Vy,2)
            wS = 1.0 - (yv[j]-ymin)/(ymax-ymin)
            Vy[i,j] = wS * BC.Vy.Dir_S[i] + (1.0-wS) * BC.Vy.Dir_N[i]
            if j>1 && j<size(Vy,2)
                wW = 1.0 - (xce[i]-xmin)/(xmax-xmin)
                Vy[i,j] += wW * BC.Vy.Dir_W[j] + (1.0-wW) * BC.Vy.Dir_E[j]
            end
        end
    end
    return
end
export SetInitialVelocity!

########

@views function SetInitialVelocity2!( Vx::Matrix{Float64}, Vy::Matrix{Float64}, BC::BoundaryConditions, dom::ModelDomain )

    ncx, ncy = dom.ncx, dom.ncy

    # Initial guess - Vx
    Vx[1    ,:] .= BC.Vx.Dir_W
    Vx[end  ,:] .= BC.Vx.Dir_E
    Vx[:,    2] .= BC.Vx.Dir_S
    Vx[:,end-1] .= BC.Vx.Dir_N
    Vy[2    ,:] .= BC.Vy.Dir_W
    Vy[end-1,:] .= BC.Vy.Dir_E
    Vy[:,    1] .= BC.Vy.Dir_S
    Vy[:,  end] .= BC.Vy.Dir_N

    Vx_ini = zeros(ncx+1,ncy+2)
    Vx_ini .= Vx
    Vy_ini = zeros(ncx+2,ncy+1)
    Vy_ini .= Vy

    for i=1:size(Vx,1)
        for j=1:size(Vx,2)
            wW = 1.0 - (dom.xv[i]-dom.xmin)/(dom.xmax-dom.xmin)
            Vx[i,j] = wW * BC.Vx.Dir_W[j] + (1.0-wW) * BC.Vx.Dir_E[j]
            if i>1 && i<size(Vx,1)
                wS = 1.0 - (dom.yce[j]-dom.ymin)/(dom.ymax-dom.ymin)
                Vx[i,j] += wS * BC.Vx.Dir_S[i] + (1.0-wS) * BC.Vx.Dir_N[i]
            end
        end
    end

    # Initial guess - Vy
    for i=1:size(Vy,1)
        for j=1:size(Vy,2)
            wS = 1.0 - (dom.yv[j]-dom.ymin)/(dom.ymax-dom.ymin)
            Vy[i,j] = wS * BC.Vy.Dir_S[i] + (1.0-wS) * BC.Vy.Dir_N[i]
            if j>1 && j<size(Vy,2)
                wW = 1.0 - (dom.xce[i]-dom.xmin)/(dom.xmax-dom.xmin)
                Vy[i,j] += wW * BC.Vy.Dir_W[j] + (1.0-wW) * BC.Vy.Dir_E[j]
            end
        end
    end
    return
end
export SetInitialVelocity2!