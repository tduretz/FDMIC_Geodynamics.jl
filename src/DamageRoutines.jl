##############

@views function DamageAssembly( f::Fields2D, Th_BC::Thermal_BC, dom::ModelDomain, mat::MaterialParameters, dt::Float64 )
    ncx, ncy, dx, dy = dom.ncx, dom.ncy, dom.dx, dom.dy
    Num       = reshape(1:ncx*ncy,ncx,ncy)
    iC        = Num
    iW        =  ones(size(Num)); iW[2:end-0,:] = Num[1:end-1,:]
    iE        =  ones(size(Num)); iE[1:end-1,:] = Num[2:end-0,:]
    iS        =  ones(size(Num)); iS[:,2:end-0] = Num[:,1:end-1]
    iN        =  ones(size(Num)); iN[:,1:end-1] = Num[:,2:end-0]
    kW        = ones(size(Num));#f.GDam .* f.lDam
    kE        = ones(size(Num));#f.GDam .* f.lDam
    kS        = ones(size(Num));#f.GDam .* f.lDam
    kN        = ones(size(Num));#f.GDam .* f.lDam
    BCW       = zeros(size(Num)); BCW[1,:]   .= Th_BC.type_W   
    BCE       = zeros(size(Num)); BCE[end,:] .= Th_BC.type_E   
    BCS       = zeros(size(Num)); BCS[:,1]   .= Th_BC.type_S   
    BCN       = zeros(size(Num)); BCN[:,end] .= Th_BC.type_S  
    # Pinned
    pC        = Th_BC.Pinned
    pW        = zeros(size(Num)); pW[2:end-0,:] = pC[1:end-1,:]
    pE        = zeros(size(Num)); pE[1:end-1,:] = pC[2:end-0,:]
    pS        = zeros(size(Num)); pS[:,2:end-0] = pC[:,1:end-1]
    pN        = zeros(size(Num)); pN[:,1:end-1] = pC[:,2:end-0]
    # Periodic stencil
    if Th_BC.type_W==0 && Th_BC.type_E==0
        println("Periodic Poisson stencil (Damage)")
        iW[  1,:] = Num[end,:]
        iE[end,:] = Num[  1,:]
    end
    # Matrix coeffcients
    fW     = 1.0 .+ 1.0.*(BCW.==1) .- 1.0.*(BCW.==2)
    fE     = 1.0 .+ 1.0.*(BCE.==1) .- 1.0.*(BCE.==2)
    fS     = 1.0 .+ 1.0.*(BCS.==1) .- 1.0.*(BCS.==2)
    fN     = 1.0 .+ 1.0.*(BCN.==1) .- 1.0.*(BCN.==2)
    cC     =  fN.*kN./dy^2 .+ fS.*kS./dy^2 .+ fE.*kE./dx^2 .+ fW.*kW./dx^2 .+ 2.0*f.We ./ (mat.gDam * mat.lDam) .+  1.0 ./mat.lDam.^2 .+ mat.eDam ./ (dt.*mat.gDam .* mat.lDam)
    cW     = -kW./dx^2; cW[BCW.!=0] .= 0.0
    cE     = -kE./dx^2; cE[BCE.!=0] .= 0.0       
    cS     = -kS./dy^2; cS[BCS.!=0] .= 0.0  
    cN     = -kN./dy^2; cN[BCN.!=0] .= 0.0  
    # Pinning
    cC[pC.==1] .= 1.0; cW[pC.==1] .= 0.0; cE[pC.==1] .= 0.0; cS[pC.==1] .= 0.0; cN[pC.==1] .= 0.0
    cW[pW.==1] .= 0.0
    cE[pE.==1] .= 0.0
    cS[pS.==1] .= 0.0
    cN[pN.==1] .= 0.0
    # Sparse matrix
    I = [iC[:]; iC[:]; iC[:]; iC[:]; iC[:]]
    J = [iC[:]; iW[:]; iE[:]; iS[:]; iN[:]]
    V = [cC[:]; cW[:]; cE[:]; cS[:]; cN[:]]
    K = sparse( I, J, V)
    droptol!(K, 1e-13)
    # display(UnicodePlots.spy(K))
    return K
end
export DamageAssembly

##############

function DamageResidual!( f::Fields2D, Th_BC::Thermal_BC, dom::ModelDomain, mat::MaterialParameters, dt::Float64 )
    """ Damage residual """
    ncx, ncy, dx, dy = dom.ncx, dom.ncy, dom.dx, dom.dy
    T = f.phiDam
    # Field with ghost values
    Tex = zeros(ncx+2, ncy+2)
    Tex[2:end-1,2:end-1] .= T
    Tex[  1,2:end-1]     .= (Th_BC.type_W==0) * T[end,:] .+ (Th_BC.type_W==2) * T[  1,:] .+ (Th_BC.type_W==1) * (2.0*Th_BC.Dir_W[:] - T[  1,:])
    Tex[end,2:end-1]     .= (Th_BC.type_E==0) * T[  1,:] .+ (Th_BC.type_E==2) * T[end,:] .+ (Th_BC.type_E==1) * (2.0*Th_BC.Dir_E[:] - T[end,:])
    Tex[2:end-1,  1]     .=                                 (Th_BC.type_S==2) * T[:,  1] .+ (Th_BC.type_S==1) * (2.0*Th_BC.Dir_S[:] - T[:,  1])
    Tex[2:end-1,end]     .=                                 (Th_BC.type_N==2) * T[:,end] .+ (Th_BC.type_N==1) * (2.0*Th_BC.Dir_N[:] - T[:,end])
    # Compute residual
    for i = 1:ncx  
        for j = 1:ncy
            pinned = Th_BC.Pinned[i,j]==1
            qW = (Tex[i+1,j+1] - Tex[i+0,j+1])/dx
            qE = (Tex[i+2,j+1] - Tex[i+1,j+1])/dx
            qS = (Tex[i+1,j+1] - Tex[i+1,j+0])/dy
            qN = (Tex[i+1,j+2] - Tex[i+1,j+1])/dy
            f.FDam[i,j] =  (pinned==0) * ( f.phiDam[i,j] * ( 2*f.We[i,j] / (mat.gDam*mat.lDam) +  1.0/mat.lDam^2 ) - 2.0*f.We[i,j] / (mat.gDam*mat.lDam) + mat.eDam*(f.phiDam[i,j]-f.phiDam0[i,j]) / (dt*mat.gDam*mat.lDam) - ( (qE-qW)/dx + (qN-qS)/dy ) ) + (pinned==1) * 0.0

            # f.FDam[i,j] =  1.0/(f.GDam[i,j]*f.lDam[i,j]) * (-2f.We[i,j] + ( f.GDam[i,j]/f.lDam[i,j] + 2f.We[i,j] )*f.phiDam[i,j]) - ( (qE-qW)/dx + (qN-qS)/dy ) 
        end
    end
    # println(norm(2f.We))
    # println(norm(f.FDam))
    # error("dd")
end
export DamageResidual!

##############

@views function InitialDamageAssembly( f::Fields2D, Th_BC::Thermal_BC, dom::ModelDomain, mat::MaterialParameters )
    ncx, ncy, dx, dy = dom.ncx, dom.ncy, dom.dx, dom.dy
    Num       = reshape(1:ncx*ncy,ncx,ncy)
    iC        = Num
    iW        =  ones(size(Num)); iW[2:end-0,:] = Num[1:end-1,:]
    iE        =  ones(size(Num)); iE[1:end-1,:] = Num[2:end-0,:]
    iS        =  ones(size(Num)); iS[:,2:end-0] = Num[:,1:end-1]
    iN        =  ones(size(Num)); iN[:,1:end-1] = Num[:,2:end-0]
    kW        = ones(size(Num));
    kE        = ones(size(Num));
    kS        = ones(size(Num));
    kN        = ones(size(Num));
    BCW       = zeros(size(Num)); BCW[1,:]   .= Th_BC.type_W   
    BCE       = zeros(size(Num)); BCE[end,:] .= Th_BC.type_E   
    BCS       = zeros(size(Num)); BCS[:,1]   .= Th_BC.type_S   
    BCN       = zeros(size(Num)); BCN[:,end] .= Th_BC.type_S  
    # Pinned
    pC        = Th_BC.Pinned
    pW        = zeros(size(Num)); pW[2:end-0,:] = pC[1:end-1,:]
    pE        = zeros(size(Num)); pE[1:end-1,:] = pC[2:end-0,:]
    pS        = zeros(size(Num)); pS[:,2:end-0] = pC[:,1:end-1]
    pN        = zeros(size(Num)); pN[:,1:end-1] = pC[:,2:end-0]
    # Periodic stencil
    if Th_BC.type_W==0 && Th_BC.type_E==0
        println("Periodic Poisson stencil (Damage)")
        iW[  1,:] = Num[end,:]
        iE[end,:] = Num[  1,:]
    end
    # Matrix coeffcients
    fW     = 1.0 .+ 1.0.*(BCW.==1) .- 1.0.*(BCW.==2)
    fE     = 1.0 .+ 1.0.*(BCE.==1) .- 1.0.*(BCE.==2)
    fS     = 1.0 .+ 1.0.*(BCS.==1) .- 1.0.*(BCS.==2)
    fN     = 1.0 .+ 1.0.*(BCN.==1) .- 1.0.*(BCN.==2)
    cC     =  fN.*kN./dy^2 .+ fS.*kS./dy^2 .+ fE.*kE./dx^2 .+ fW.*kW./dx^2 .+ 1.0 / mat.lDam^2
    cW     = -kW./dx^2; cW[BCW.!=0] .= 0.0
    cE     = -kE./dx^2; cE[BCE.!=0] .= 0.0       
    cS     = -kS./dy^2; cS[BCS.!=0] .= 0.0  
    cN     = -kN./dy^2; cN[BCN.!=0] .= 0.0      
    # Pinning
    cC[pC.==1] .= 1.0; cW[pC.==1] .= 0.0; cE[pC.==1] .= 0.0; cS[pC.==1] .= 0.0; cN[pC.==1] .= 0.0
    cW[pW.==1] .= 0.0
    cE[pE.==1] .= 0.0
    cS[pS.==1] .= 0.0
    cN[pN.==1] .= 0.0
    # Sparse matrix
    I = [iC[:]; iC[:]; iC[:]; iC[:]; iC[:]]
    J = [iC[:]; iW[:]; iE[:]; iS[:]; iN[:]]
    V = [cC[:]; cW[:]; cE[:]; cS[:]; cN[:]]
    K = sparse( I, J, V)
    droptol!(K, 1e-13)
    # display(UnicodePlots.spy(K))
    return K
end
export InitialDamageAssembly

##############

function InitialDamageResidual!( f::Fields2D, Th_BC::Thermal_BC, dom::ModelDomain, mat::MaterialParameters )
    """ Damage residual """
    ncx, ncy, dx, dy = dom.ncx, dom.ncy, dom.dx, dom.dy
    T = f.phiDam
    # Field with ghost values
    Tex = zeros(ncx+2, ncy+2)
    Tex[2:end-1,2:end-1] .= T
    Tex[  1,2:end-1]     .= (Th_BC.type_W==0) * T[end,:] .+ (Th_BC.type_W==2) * T[  1,:] .+ (Th_BC.type_W==1) * (2.0*Th_BC.Dir_W[:] - T[  1,:])
    Tex[end,2:end-1]     .= (Th_BC.type_E==0) * T[  1,:] .+ (Th_BC.type_E==2) * T[end,:] .+ (Th_BC.type_E==1) * (2.0*Th_BC.Dir_E[:] - T[end,:])
    Tex[2:end-1,  1]     .=                                 (Th_BC.type_S==2) * T[:,  1] .+ (Th_BC.type_S==1) * (2.0*Th_BC.Dir_S[:] - T[:,  1])
    Tex[2:end-1,end]     .=                                 (Th_BC.type_N==2) * T[:,end] .+ (Th_BC.type_N==1) * (2.0*Th_BC.Dir_N[:] - T[:,end])
    # Compute residual
    for i = 1:ncx  
        for j = 1:ncy
            pinned = Th_BC.Pinned[i,j]==1
            qW = (Tex[i+1,j+1] - Tex[i+0,j+1])/dx
            qE = (Tex[i+2,j+1] - Tex[i+1,j+1])/dx
            qS = (Tex[i+1,j+1] - Tex[i+1,j+0])/dy
            qN = (Tex[i+1,j+2] - Tex[i+1,j+1])/dy
            if pinned==1
                f.FDam[i,j] = 0.0
            else
                f.FDam[i,j] = f.phiDam[i,j] / mat.lDam^2 - ( (qE-qW)/dx + (qN-qS)/dy ) 
                # f.FDam[i,j] = f.phiDam[i,j] / mat.lDam[i,j]^2 - ( (qE-qW)/dx + (qN-qS)/dy ) 
            end
        end
    end
end
export InitialDamageResidual!