

@views function StrainRate!(div,Eps,Vx,Vy,dx,dy)
    div    .= diff(Vx[:,2:end-1],dims=1)./dx + diff(Vy[2:end-1,:],dims=2)./dy 
    Eps.xx .= diff(Vx[:,2:end-1],dims=1)./dx .- 1.0/3.0 .* div
    Eps.yy .= diff(Vy[2:end-1,:],dims=2)./dy .- 1.0/3.0 .* div
    Eps.zz .= -(Eps.xx + Eps.yy)
    Eps.xy .= 0.5.*( diff(Vx,dims=2)./dy .+ diff(Vy,dims=1)/dx ) 
    Exy     = 0.25*(Eps.xy[1:end-1,1:end-1] .+ Eps.xy[1:end-1,2:end-0] .+ Eps.xy[2:end-0,1:end-1] .+ Eps.xy[2:end-0,2:end-0])
    Eps.II .= sqrt.(0.5*(Eps.xx.^2 .+ Eps.xx.^2 .+ Eps.zz.^2) .+ Exy.^2)
end

########

@views function Stress!(Tau,Eps,etac,etav)
    Tau.xx .= 2.0.*etac.*Eps.xx 
    Tau.yy .= 2.0.*etac.*Eps.yy
    Tau.zz .= 2.0.*etac.*Eps.zz  
    Tau.xy .= 2.0.*etav.*Eps.xy 
    Txy     = 0.25*(Tau.xy[1:end-1,1:end-1] .+ Tau.xy[1:end-1,2:end-0] .+ Tau.xy[2:end-0,1:end-1] .+ Tau.xy[2:end-0,2:end-0])
    Tau.II .= sqrt.(0.5*(Tau.xx.^2 .+ Tau.xx.^2 .+ Tau.zz.^2) .+ Txy.^2)
end

########

@views function Residuals!(Fx,Fy,Fp,Tau,Pc,div,BC,ncx,ncy,dx,dy)
    Fx .= 0.0
    Fy .= 0.0
    Fp .= 0.0
    if BC.periodix==0
        Fx[2:end-1,:] .= -( diff(Tau.xx .- Pc ,dims=1)./dx .+ diff(Tau.xy[2:end-1,:],dims=2)/dy )
    else
        Sxx_ex             = zeros(ncx+1,ncy)
        Sxx_ex[2:end-0,:] .= -Pc .+ Tau.xx
        Sxx_ex[      1,:] .= -Pc[end,:] .+ Tau.xx[end,:]
        #Sxx_ex[    end,:] .= -Pc[  1,:] .+ Tau.xx[  1,:] # Do not assemble last column
        Fx .= -( diff(Sxx_ex ,dims=1)./dx .+ diff(Tau.xy[1:end-1,:],dims=2)/dy )
    end
    Fy[:,2:end-1] .= -( diff(Tau.yy .- Pc ,dims=2)./dy .+ diff(Tau.xy[:,2:end-1],dims=1)/dx )
    Fp            .= -div
    # # For periodic
    # if BC.periodix==1
    #     Fx = Fx[1:end-1,:]
    # end
end

########

@views function SetBCs( Vx, Vy, BC )
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

########

function NumberingStokes(BC, ncx, ncy)
    # Numbering
    if BC.periodix==0
        NumVx     = collect(reshape(1:(ncx+1)*ncy,ncx+1,ncy))
    else
        NumVx             = zeros(Int64,(ncx+1),ncy)
        NumVx[1:end-1,:] .= collect(reshape(1:(ncx+0)*ncy,ncx,ncy))
        NumVx[end,:]     .= NumVx[1,:]
    end
    NumVy     = collect(reshape(1:ncx*(ncy+1),ncx,ncy+1) .+ maximum(NumVx))
    NumP      = collect(reshape(1:(ncx)*ncy,ncx,ncy))
    return NumVx, NumVy, NumP
end

########

@views function StokesAssembly( BC, NumVx, NumVy, NumP, etac, etav, DirScale, dx, dy )

    # Connectivity
    iVxC      = NumVx
    iVxW      =  ones(Int64, size(NumVx)); iVxW[2:end-0,: ] = NumVx[1:end-1,:]
    iVxE      =  ones(Int64, size(NumVx)); iVxE[1:end-1,: ] = NumVx[2:end-0,:]        
    iVxS      =  ones(Int64, size(NumVx)); iVxS[: ,2:end-0] = NumVx[:,1:end-1]
    iVxN      =  ones(Int64, size(NumVx)); iVxN[: ,1:end-1] = NumVx[:,2:end-0]
    iVySW     =  ones(Int64, size(NumVx)); iVySW[2:end-0,:] = NumVy[:,1:end-1]
    iVySE     =  ones(Int64, size(NumVx)); iVySE[1:end-1,:] = NumVy[:,1:end-1]
    iVyNW     =  ones(Int64, size(NumVx)); iVyNW[2:end-0,:] = NumVy[:,2:end-0]
    iVyNE     =  ones(Int64, size(NumVx)); iVyNE[1:end-1,:] = NumVy[:,2:end-0]
    iPW       =  ones(Int64, size(NumVx)); iPW[2:end-0,:]   = NumP
    iPE       =  ones(Int64, size(NumVx)); iPE[1:end-1,:]   = NumP

    # Viscosity coefficients
    etaW      = zeros(size(NumVx)); etaW[2:end-0,:] = etac[1:end-0,:]
    etaE      = zeros(size(NumVx)); etaE[1:end-1,:] = etac[1:end-0,:]
    etaS      = zeros(size(NumVx)); etaS[:,1:end-0] = etav[:,1:end-1] 
    etaN      = zeros(size(NumVx)); etaN[:,1:end-0] = etav[:,2:end-0]

    if BC.periodix==1
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
    cPW   = -1.0/dx .*  ones(size(NumVx))
    cPE   =  1.0/dx .*  ones(size(NumVx))

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
    cPS   = -1.0/dy .* ones(size(NumVy)); cPS[:,  1] .= 0.0;  cPS[:,end] .= 0.0
    cPN   =  1.0/dy .* ones(size(NumVy)); cPN[:,  1] .= 0.0;  cPN[:,end] .= 0.0

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
        iVxC  = iVxC[1:end-1,:];  cVxC  = cVxC[1:end-1,:]
        iVxW  = iVxW[1:end-1,:];  cVxW  = cVxW[1:end-1,:]
        iVxE  = iVxE[1:end-1,:];  cVxE  = cVxE[1:end-1,:]
        iVxS  = iVxS[1:end-1,:];  cVxS  = cVxS[1:end-1,:]
        iVxN  = iVxN[1:end-1,:];  cVxN  = cVxN[1:end-1,:]
        iVySW = iVySW[1:end-1,:]; cVySW = cVySW[1:end-1,:]
        iVySE = iVySE[1:end-1,:]; cVySE = cVySE[1:end-1,:]
        iVyNW = iVyNW[1:end-1,:]; cVyNW = cVyNW[1:end-1,:]
        iVyNE = iVyNE[1:end-1,:]; cVyNE = cVyNE[1:end-1,:]
        iPW   = iPW[1:end-1,:];   cPW   = cPW[1:end-1,:]
        iPE   = iPE[1:end-1,:];   cPE   = cPE[1:end-1,:]
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
    # Sparse matrix Kpu
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

########

@views function StokesSolver!(Vx,Vy,Pc,NumVx,NumVy,NumP,fu,fp,Kuu,Kup,Kpu,etac,gamma,solver)

    if solver == 0
        # Slightly compressible pressure bloack
        cP   = 1.0./gamma.*ones(Float64,size(etac,1),size(etac,2))
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
    else
        DecoupledSolver!(Vx,Vy,Pc,NumVx,NumVy,NumP,fu,fp,Kuu,Kup,Kpu,etac,gamma)
    end
    return
end

########

@views function DecoupledSolver!(Vx,Vy,Pc,NumVx,NumVy,NumP,fu,fp,Kuu,Kup,Kpu,etac,gamma)
    # Decoupled solve
    ndofu = size(Kup,1)
    ndofp = size(Kup,2)
    coef  = gamma*ones(length(etac))#.*etac[:]
    Kppi  = spdiagm(coef)
    Kuusc = Kuu - Kup*(Kppi*Kpu) # OK
    PC    =  0.5*(Kuusc + Kuusc') 
    t = @elapsed Kf    = cholesky(Hermitian(PC),check = false)
    @printf("Cholesky took = %02.2e s\n", t)
    u     = zeros(ndofu, 1)
    ru    = zeros(ndofu, 1)
    fusc  = zeros(ndofu, 1)
    p     = zeros(ndofp, 1)
    rp    = zeros(ndofp, 1)
    # Iterations
    for rit=1:10
        ru   .= fu - Kuu*u - Kup*p;
        rp   .= fp - Kpu*u;
        @printf("  --> Powell-Hestenes Iteration %02d\n  Momentum res.   = %2.2e\n  Continuity res. = %2.2e\n", rit, norm(ru)/sqrt(length(ru)), norm(rp)/sqrt(length(rp)))
        if norm(ru)/(length(ru)) < 1e-10 && norm(rp)/(length(ru)) < 1e-10
            break
        end
        fusc .=  fu  - Kup*(Kppi*fp + p)
        u    .= Kf\fusc
        p   .+= Kppi*(fp - Kpu*u)
    end
    Vx[:,2:end-1] .+= u[NumVx]
    Vy[2:end-1,:] .+= u[NumVy]
    Pc            .+= p[NumP]
end

########

@views function SetInitialVelocity!( Vx, Vy, BC, xv, yv, xmin, xmax, ymin, ymax, ncx, ncy )

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
            # if i>1 & i<size(Vx,1)
            #     wS = 1.0 - (yce[j]-ymin)/(ymax-ymin)
            #     Vx[i,j] += wS * BC.Vx.Dir_S[i] + (1.0-wS) * BC.Vx.Dir_N[i]
            # end
        end
    end

    # Initial guess - Vy
    for i=1:size(Vy,1)
        for j=1:size(Vy,2)
            wS = 1.0 - (yv[j]-ymin)/(ymax-ymin)
            Vy[i,j] = wS * BC.Vy.Dir_S[i] + (1.0-wS) * BC.Vy.Dir_N[i]
            # if j>1 & j<size(Vy,2)
            #     wW = 1.0 - (xce[i]-xmin)/(xmax-xmin)
            #     Vy[i,j] += wW * BC.Vy.Dir_W[j] + (1.0-wW) * BC.Vy.Dir_E[j]
            # end
        end
    end
    return
end