function VerticesToCentroids!( kc::Matrix{Float64}, kv::Matrix{Float64} )
    @tturbo for j=1:size(kc,2)
        for i=1:size(kc,1)
            kc[i,j] = 0.25*(kv[i,j] + kv[i+1,j] + kv[i,j+1] + kv[i+1,j+1])
        end
    end
end
export VerticesToCentroids!

##############

function CentroidsToVertices!( kv::Matrix{Float64}, kc::Matrix{Float64}, ncx::Int64, ncy::Int64, BC::BoundaryConditions )

    @tturbo kv .= 0.0

    # First make extended arry using BC's
    kc_ex   = zeros( ncx+2, ncy+2)
    @tturbo kc_ex[2:end-1,2:end-1] .= kc
    if BC.periodix==1 # Periodic case
        kc_ex[  1,2:end-1] .= kc[end,:] 
        kc_ex[end,2:end-1] .= kc[  1,:]
    else # other cases
        kc_ex[  1,2:end-1] .= kc[  1,:] 
        kc_ex[end,2:end-1] .= kc[end,:]
    end
    kc_ex[  :,  1] .= kc_ex[  :,    2]
    kc_ex[  :,end] .= kc_ex[  :,end-1]

    # Conductivity averaging
    VerticesToCentroids!( kv, kc_ex )
end
export CentroidsToVertices!

##############

function AverageConductivity( kc::Matrix{Float64}, ncx::Int64, ncy::Int64, Th_BC::Thermal_BC )

    # Security
    if Th_BC.type_W==0 && Th_BC.type_E!=0
        error("Both sides should be periodic!")
    end
    kx        = zeros( ncx+1, ncy  )
    ky        = zeros( ncx  , ncy+1)
    kc_ex     = zeros( ncx+2, ncy  )
    kc_ex[2:end-1,:] .= kc
    if Th_BC.type_W==0 && Th_BC.type_E==0 # Periodic case
        kc_ex[  1,:] .= kc[end,:] 
        kc_ex[end,:] .= kc[  1,:]
    else # other cases
        kc_ex[  1,:] .= kc[  1,:] 
        kc_ex[end,:] .= kc[end,:]
    end

    # Conductivity averaging x
    @tturbo for j=1:ncy 
        for i=1:ncx+1
            kx[i,j] += 0.5*(kc_ex[i,j] + kc_ex[i+1,j])
        end
    end

    kc_ex     = zeros( ncx, ncy+2)
    kc_ex[:,2:end-1] .= kc
    kc_ex[:,  1] .= kc[:,  1]
    kc_ex[:,end] .= kc[:,end]
    
    # Conductivity averaging y
    @tturbo for i=1:ncx
        for j=1:ncy+1
            ky[i,j] += 0.5*(kc_ex[i,j] + kc_ex[i,j+1])
        end
    end
    
    return kx, ky
end
export AverageConductivity

##############

@views function GenerateMesh( xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64, ncx::Int64, ncy::Int64 )
    dx, dy    = (xmax-xmin)/ncx, (ymax-ymin)/ncy
    xc        = LinRange( xmin+dx/2, xmax-dx/2, ncx  )
    yc        = LinRange( ymin+dy/2, ymax-dy/2, ncy  )
    xce       = LinRange( xmin-dx/2, xmax+dx/2, ncx+2)
    yce       = LinRange( ymin-dy/2, ymax+dy/2, ncy+2)
    xv        = LinRange( xmin     , xmax     , ncx+1)
    yv        = LinRange( ymin     , ymax     , ncy+1)
    return dx, dy, xc, yc, xce, yce, xv, yv
end
export GenerateMesh

##############

@views function GenerateMesh!( dom::ModelDomain )
    """
    This functions generate x and y coordinates of mesh centroids and vertices, it also computes the grid spacing
    """
    dom.dx, dom.dy    = (dom.xmax-dom.xmin)/dom.ncx, (dom.ymax-dom.ymin)/dom.ncy
    dom.xc        = LinRange( dom.xmin+dom.dx/2, dom.xmax-dom.dx/2, dom.ncx  )
    dom.yc        = LinRange( dom.ymin+dom.dy/2, dom.ymax-dom.dy/2, dom.ncy  )
    dom.xce       = LinRange( dom.xmin-dom.dx/2, dom.xmax+dom.dx/2, dom.ncx+2)
    dom.yce       = LinRange( dom.ymin-dom.dy/2, dom.ymax+dom.dy/2, dom.ncy+2)
    dom.xv        = LinRange( dom.xmin     , dom.xmax     , dom.ncx+1)
    dom.yv        = LinRange( dom.ymin     , dom.ymax     , dom.ncy+1)
    return 
end
export GenerateMesh!

##############

@views function PureShearBoxUpdate!( dom::ModelDomain, BC::BoundaryConditions, params::ModelParameters )
    dom.xmin    += mean(BC.Vx.Dir_W) * params.dt
    dom.xmax    += mean(BC.Vx.Dir_E) * params.dt
    dom.ymin    += mean(BC.Vy.Dir_S) * params.dt
    dom.ymax    += mean(BC.Vy.Dir_N) * params.dt            
    BC.Vx.Dir_W .=-params.Ebg*dom.xmin*ones(dom.ncy+2)
    BC.Vx.Dir_E .=-params.Ebg*dom.xmax*ones(dom.ncy+2)            
    BC.Vy.Dir_S .= params.Ebg*dom.ymin*ones(dom.ncx+2)
    BC.Vy.Dir_N .= params.Ebg*dom.ymax*ones(dom.ncx+2)
    dom.L = dom.xmax - dom.xmin
    println("Box pure shear deformation: ", (dom.L-dom.L0)/dom.L0*100, "%" )
    GenerateMesh!( dom )
    return dom.dx, dom.dy
end
export PureShearBoxUpdate!