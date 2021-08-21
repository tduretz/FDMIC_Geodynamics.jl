function VerticesToCentroids!( kc::Matrix{Float64}, kv::Matrix{Float64} )
    @tturbo for j=1:size(kc,2)
        for i=1:size(kc,1)
            kc[i,j] += 0.25*(kv[i,j] + kv[i+1,j] + kv[i,j+1] + kv[i+1,j+1])
        end
    end
end


function CentroidsToVertices!( kv, kc, ncx, ncy, BC )

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

##############

function AverageConductivity( kc, ncx, ncy, Th_BC )

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

##############

@views function GenerateMesh( xmin, xmax, ymin, ymax, ncx, ncy )
    dx, dy    = (xmax-xmin)/ncx, (ymax-ymin)/ncy
    xc        = LinRange( xmin+dx/2, xmax-dx/2, ncx  )
    yc        = LinRange( ymin+dy/2, ymax-dy/2, ncy  )
    xce       = LinRange( xmin-dx/2, xmax+dx/2, ncx+2)
    yce       = LinRange( ymin-dy/2, ymax+dy/2, ncy+2)
    xv        = LinRange( xmin     , xmax     , ncx+1)
    yv        = LinRange( ymin     , ymax     , ncy+1)
    return dx, dy, xc, yc, xce, yce, xv, yv
end

##############

@views function PureShearBoxUpdate!( xmin, xmax, ymin, ymax, BC, Ebg, dt, L0, ncx, ncy )
    

    # L = xmax - xmin
    # println("Box pure shear deformation: ", (L-L0)/L0*100, "%" )
    return
end