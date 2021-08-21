##############

@views function ThermalAssembly( Cp, rho, kx, ky, dx, dy, dt, steady, Th_BC, ncx, ncy )

    Num       = reshape(1:ncx*ncy,ncx,ncy)
    iC        = Num
    iW        =  ones(size(Num)); iW[2:end-0,:] = Num[1:end-1,:]
    iE        =  ones(size(Num)); iE[1:end-1,:] = Num[2:end-0,:]
    iS        =  ones(size(Num)); iS[:,2:end-0] = Num[:,1:end-1]
    iN        =  ones(size(Num)); iN[:,1:end-1] = Num[:,2:end-0]
    kW        = zeros(size(Num)); kW[:,:] = kx[1:end-1,:]
    kE        = zeros(size(Num)); kE[:,:] = kx[2:end-0,:]
    kS        = zeros(size(Num)); kS[:,:] = ky[:,1:end-1] 
    kN        = zeros(size(Num)); kN[:,:] = ky[:,2:end-0]
    BCW       = zeros(size(Num)); BCW[1,:]   .= Th_BC.type_W   
    BCE       = zeros(size(Num)); BCE[end,:] .= Th_BC.type_E   
    BCS       = zeros(size(Num)); BCS[:,1]   .= Th_BC.type_S   
    BCN       = zeros(size(Num)); BCN[:,end] .= Th_BC.type_S  
    # Periodic stencil
    if Th_BC.type_W==0 && Th_BC.type_E==0
        iW[  1,:] = Num[end,:]
        iE[end,:] = Num[  1,:]
    end
    # Matrix coeffcients
    fW     = 1.0 .+ 1.0.*(BCW.==1) .- 1.0.*(BCW.==2)
    fE     = 1.0 .+ 1.0.*(BCE.==1) .- 1.0.*(BCE.==2)
    fS     = 1.0 .+ 1.0.*(BCS.==1) .- 1.0.*(BCS.==2)
    fN     = 1.0 .+ 1.0.*(BCN.==1) .- 1.0.*(BCN.==2)
    cC     =  (steady!=1) .* Cp.*rho./dt .+ fN.*kN./dy^2 .+ fS.*kS./dy^2 .+ fE.*kE./dx^2 .+ fW.*kW./dx^2
    cW     = -kW./dx^2; cW[BCW.!=0] .= 0.0
    cE     = -kE./dx^2; cE[BCE.!=0] .= 0.0       
    cS     = -kS./dx^2; cS[BCS.!=0] .= 0.0  
    cN     = -kN./dx^2; cN[BCN.!=0] .= 0.0  
    
    # # Matrix coeffcients
    # cC        = zeros(size(Num))
    # cW        = zeros(size(Num))
    # cE        = zeros(size(Num))
    # cS        = zeros(size(Num))
    # cN        = zeros(size(Num))
    # @threads for i = 1:ncx # vectorisation does not work
    #     for j = 1:ncy
    #         # Setup BCs
    #         DirW    = (i==1) & (Th_BC.type_W==1)
    #         NeuW    = (i==1) & (Th_BC.type_W==2)
    #         fW      = 1.0 + 1.0*DirW - 1.0*NeuW 
    #         DirE    = (i==ncx) & (Th_BC.type_E==1)
    #         NeuE    = (i==ncx) & (Th_BC.type_E==2)
    #         fE      = 1.0 + 1.0*DirE - 1.0*NeuE 
    #         DirS    = (j==1) & (Th_BC.type_S==1)
    #         NeuS    = (j==1) & (Th_BC.type_S==2)
    #         fS      = 1.0 + 1.0*DirS - 1.0*NeuS 
    #         DirN    = (j==ncy) & (Th_BC.type_N==1)
    #         NeuN    = (j==ncy) & (Th_BC.type_N==2)
    #         fN      = 1.0 + 1.0*DirN - 1.0*NeuN
    #         # Coefficients
    #         cC[i,j] =  (steady!=1) * Cp[i,j]*rho[i,j]/dt + fN*kN[i,j]/dy^2 + fS*kS[i,j]/dy^2 + fE*kE[i,j]/dx^2 + fW*kW[i,j]/dx^2
    #         cW[i,j] = -(i>1  ) * kW[i,j]/dx^2 
    #         cE[i,j] = -(i<ncx) * kE[i,j]/dx^2 
    #         cS[i,j] = -(j>  1) * kS[i,j]/dy^2 
    #         cN[i,j] = -(j<ncy) * kN[i,j]/dy^2 
    #     end
    # end
    
    # Sparse matrix
    I = [iC[:]; iC[:]; iC[:]; iC[:]; iC[:]]
    J = [iC[:]; iW[:]; iE[:]; iS[:]; iN[:]]
    V = [cC[:]; cW[:]; cE[:]; cS[:]; cN[:]]
    K = sparse( I, J, V)
    droptol!(K, 1e-8)
    # display(UnicodePlots.spy(K))
    return K
end

##############

@views function ThermalSolve!( T, K, F, ncx, ncy )
    # Solve
    b     = zeros(ncx*ncy,1); b[:,1] = F[:]
    Kchol = cholesky(K)
    dT    = Kchol\b
    T   .-= reshape(dT,ncx,ncy)
end

##############

function ThermalResidual(T, T0, rho, Cp, H, kx, ky, dx, dy, dt, steady, Th_BC, ncx, ncy)
    """ Thermal residuals """
    Tex = zeros(ncx+2, ncy+2)
    Tex[2:end-1,2:end-1] .= T
    Tex[  1,2:end-1]     .= (Th_BC.type_W==0) * T[end,:] .+ (Th_BC.type_W==2) * T[  1,:] .+ (Th_BC.type_W==1) * (2.0*Th_BC.Dir_W[:] - T[  1,:])
    Tex[end,2:end-1]     .= (Th_BC.type_E==0) * T[  1,:] .+ (Th_BC.type_E==2) * T[end,:] .+ (Th_BC.type_E==1) * (2.0*Th_BC.Dir_E[:] - T[end,:])
    Tex[2:end-1,  1]     .=                                 (Th_BC.type_S==2) * T[:,  1] .+ (Th_BC.type_S==1) * (2.0*Th_BC.Dir_S[:] - T[:,  1])
    Tex[2:end-1,end]     .=                                 (Th_BC.type_N==2) * T[:,end] .+ (Th_BC.type_N==1) * (2.0*Th_BC.Dir_N[:] - T[:,end])

    # Compute residual
    F = zeros( ncx, ncy )
    @tturbo for i = 1:ncx
        for j = 1:ncy
            qW = -kx[i  ,j]*(Tex[i+1,j+1] - Tex[i+0,j+1])/dx
            qE = -kx[i+1,j]*(Tex[i+2,j+1] - Tex[i+1,j+1])/dx
            qS = -ky[i,j  ]*(Tex[i+1,j+1] - Tex[i+1,j+0])/dy
            qN = -ky[i,j+1]*(Tex[i+1,j+2] - Tex[i+1,j+1])/dy
            F[i,j] = (steady!=1) * rho[i,j]*Cp[i,j] * (T[i,j] - T0[i,j])/dt + (qE-qW)/dx + (qN-qS)/dy + H[i,j]
        end
    end
return F
end

##############