using Base.Threads

function InitialiseMarkerPositions!( p::Markers, domain::ModelDomain )
    p.nmark     = domain.ncx*domain.ncy*p.nmx*p.nmy; # total initial number of marker in grid
    p.nmark_ini = p.nmark; # total initial number of marker in grid
    p.nmark_max = 3*p.nmark; # total initial number of marker in grid
    dxm, dym    = domain.dx/p.nmx, domain.dy/p.nmy 
    xm1d        = LinRange(domain.xmin+dxm/2, domain.xmax-dxm/2, domain.ncx*p.nmx)
    ym1d        = LinRange(domain.ymin+dym/2, domain.ymax-dym/2, domain.ncy*p.nmy)
    (xmi,ymi)   = ([x for x=xm1d,y=ym1d], [y for x=xm1d,y=ym1d])
    p.x         = zeros(Float64, p.nmark_max)
    p.y         = zeros(Float64, p.nmark_max)
    p.x[1:p.nmark] .= vec(xmi)
    p.y[1:p.nmark] .= vec(ymi)
    return
end
export InitialiseMarkerPositions!

@views function CopyMarkers!(p::Markers,ineigh::Int64)
    p.phase[p.nmark] = p.phase[ineigh]
    p.T[p.nmark]     = p.T[ineigh]
end
export CopyMarkers!

@views function FindClosestNeighbour( p::Markers, mpc::Matrix{Float64}, mlist::Matrix{Array{Int64, N} where N}, new_x::Float64, new_y::Float64, dxy::Float64, i::Int64, j::Int64, ncx::Int64, ncy::Int64 )
    ineigh = 0
    dst = dxy
    # Current cell
    for im=1:Int64(mpc[i,j])
        xm = p.x[mlist[i,j][im]] 
        ym = p.y[mlist[i,j][im]] 
        d = (new_x-xm)^2 + (new_y-ym)^2
        if d<dst
            dst    = d
            ineigh = mlist[i,j][im]
        end
    end
    if ineigh==0
        # Direct neighbours: W
        if i>1 # could make a periodic case here...
            for im=1:Int64(mpc[i-1,j])
                xm = p.x[mlist[i-1,j][im]] 
                ym = p.y[mlist[i-1,j][im]] 
                d = (new_x-xm)^2 + (new_y-ym)^2
                if d<dst
                    dst    = d
                    ineigh = mlist[i-1,j][im]
                end
            end
        end
        # Direct neighbours: E
        if i<ncx # could make a periodic case here...
            for im=1:Int64(mpc[i+1,j])
                xm = p.x[mlist[i+1,j][im]] 
                ym = p.y[mlist[i+1,j][im]] 
                d = (new_x-xm)^2 + (new_y-ym)^2
                if d<dst
                    dst    = d
                    ineigh = mlist[i+1,j][im]
                end
            end
        end
        # Direct neighbours: S
        if j>1 # could make a periodic case here...
            for im=1:Int64(mpc[i,j-1])
                xm = p.x[mlist[i,j-1][im]] 
                ym = p.y[mlist[i,j-1][im]] 
                d = (new_x-xm)^2 + (new_y-ym)^2
                if d<dst
                    dst    = d
                    ineigh = mlist[i,j-1][im]
                end
            end
        end
        # Direct neighbours: S
        if j<ncy # could make a periodic case here...
            for im=1:Int64(mpc[i,j+1])
                xm = p.x[mlist[i,j+1][im]] 
                ym = p.y[mlist[i,j+1][im]] 
                d = (new_x-xm)^2 + (new_y-ym)^2
                if d<dst
                    dst    = d
                    ineigh = mlist[i,j+1][im]
                end
            end
        end
    end
    # Last check, throw error
    if (ineigh==0) error("Some cell ran empty of markers...") end
    return ineigh
end
export FindClosestNeighbour

@views function ReseedMarkers!( p::Markers, mpc::Matrix{Float64}, xc::LinRange{Float64}, yc::LinRange{Float64}, dx::Float64, dy::Float64, ncx::Int64, ncy::Int64)

    mlist       = Array{Array{Int64}}(undef,ncx,ncy)
    mlist_th    = Array{Array{Int64}}(undef,nthreads(),ncx,ncy)

    for i=1:ncx
        for j=1:ncy
            mlist[i,j] = zeros(Int64(mpc[i,j]))
        end
    end

    count = zeros(Int64,ncx,ncy)

    for k=1:p.nmark # @avx ne marche pas ici
        if p.phase[k]>=0
            # Get the column:
            dstx = p.x[k] - xc[1];
            i    = Int64(round(ceil( (dstx/dx) + 0.5)));
            # Get the line:
            dsty = p.y[k] - yc[1];
            j    = Int64(round(ceil( (dsty/dy) + 0.5)));
            # Increment cell count
            count[i,j] += 1
            mlist[i,j][count[i,j]] = k
        end
    end

    # display(mlist)
    dxy = 10*(dx*dx + dy*dy)
    for i=1:ncx
        for j=1:ncy

            if mpc[i,j] < 4

                new_x = xc[i]-dx/4.0
                new_y = yc[j]-dy/4.0
                ineigh = FindClosestNeighbour( p, mpc, mlist, new_x, new_y, dxy, i, j, ncx, ncy )
                if ineigh>0
                    p.nmark+=1
                    p.x[p.nmark] = new_x
                    p.y[p.nmark] = new_y
                    CopyMarkers!(p,ineigh)
                end

                new_x = xc[i]+dx/4.0
                new_y = yc[j]-dy/4.0
                ineigh = FindClosestNeighbour( p, mpc, mlist, new_x, new_y, dxy, i, j, ncx, ncy )
                if ineigh>0
                    p.nmark+=1
                    p.x[p.nmark] = new_x
                    p.y[p.nmark] = new_y
                    CopyMarkers!(p,ineigh)
                end

                new_x = xc[i]+dx/4.0
                new_y = yc[j]+dy/4.0
                ineigh = FindClosestNeighbour( p, mpc, mlist, new_x, new_y, dxy, i, j, ncx, ncy )
                if ineigh>0
                    p.nmark+=1
                    p.x[p.nmark] = new_x
                    p.y[p.nmark] = new_y
                    CopyMarkers!(p,ineigh)
                end

                new_x = xc[i]-dx/4.0
                new_y = yc[j]+dy/4.0
                ineigh = FindClosestNeighbour( p, mpc, mlist, new_x, new_y, dxy, i, j, ncx, ncy )
                if ineigh>0
                    p.nmark+=1
                    p.x[p.nmark] = new_x
                    p.y[p.nmark] = new_y
                    CopyMarkers!(p,ineigh)
                end

            end
        end
    end
end
export ReseedMarkers!

@views function ListMarkers( p::Markers, ncx, ncy )
liste = hcat([[Int[] for i in 1:ncx] for j in 1:ncy]...)
for k in 1:p.nmark
    if p.phase[k]>=0
        i = p.cellx[k]
        j = p.celly[k]
        push!(liste[i,j], k)
    end
end
return liste
end
export ListMarkers


function RungeKutta!(p::Markers, nmark::Int64, rkv::Matrix{Float64}, rkw::Matrix{Float64}, BC::BoundaryConditions, dt::Float64, Vx::Matrix{Float64}, Vy::Matrix{Float64}, xv::LinRange{Float64}, yv::LinRange{Float64}, xce::LinRange{Float64}, yce::LinRange{Float64}, dx::Float64, dy::Float64, ncx::Int64, ncy::Int64)
# Marker advection with 4th order Roger-Kutta
@threads for k=1:nmark
    in = p.phase[k]>=0
        x0 = p.x[k];
        y0 = p.y[k];
        vx = 0.0
        vy = 0.0
        # Roger-Gunther loop
        for rk=1:4
            # println("step", rk)
            # Interp velocity from grid
            vxm = VxFromVxNodes(Vx, k, p, xv, yce, dx, dy, ncx, ncy, 1)
            vym = VyFromVxNodes(Vy, k, p, xce, yv, dx, dy, ncx, ncy, 1)
            # Temporary RK advection steps
            p.x[k] = x0 + rkv[rk]*dt*vxm
            p.y[k] = y0 + rkv[rk]*dt*vym
            if BC.periodix==1
                if p.x[k]<xv[1]
                    # println("periW 1 ", p.x[k], vxm)
                    p.x[k]+=(xv[end]-xv[1]) 
                    # println("periW 2 ", p.x[k], vxm)
                end
                if p.x[k]>xv[end] 
                    # println("periE")
                    p.x[k]-=(xv[end]-xv[1]) 
                end
            end
            # Average final velocity 
            vx    += rkw[rk]*vxm
            vy    += rkw[rk]*vym
        end
        # Advect points
        p.x[k] = x0 + (in==1) * rkv[4]*dt*vx
        p.y[k] = y0 + (in==1) * rkv[4]*dt*vy
        if BC.periodix==1
            if p.x[k]<xv[1]
                # println("periW 1 ", p.x[k], vxm)
                p.x[k]+=(xv[end]-xv[1]) 
                # println("periW 2 ", p.x[k], vxm)
            end
            if p.x[k]>xv[end] 
                # println("periE")
                p.x[k]-=(xv[end]-xv[1]) 
            end
        end
end
end
export RungeKutta!

# Interpolation from Vx nodes to particles
@views function VxFromVxNodes(Vx::Matrix{Float64}, k::Int64, p::Markers, xv::LinRange{Float64}, yce::LinRange{Float64}, dx::Float64, dy::Float64, ncx::Int64, ncy::Int64, new::Int64)
    # Interpolate vx
    i = Int64(round(trunc( (p.x[k] -  xv[1])/dx ) + 1));
    j = Int64(round(trunc( (p.y[k] - yce[1])/dy ) + 1));
    if i<1
        i = 1
    elseif i>ncx
        i = ncx
    end
    if j<1
        j = 1;
    elseif j> ncy+1
        j = ncy+1
    end
    # Compute distances
    dxmj = p.x[k] -  xv[i]
    dymi = p.y[k] - yce[j]
    # Compute vx velocity for the top and bottom of the cell
    vxm13 = Vx[i,j  ] * (1-dxmj/dx) + Vx[i+1,j  ]*dxmj/dx
    vxm24 = Vx[i,j+1] * (1-dxmj/dx) + Vx[i+1,j+1]*dxmj/dx
    if new==1 
        if dxmj/dx>=0.5
            if i<ncx
                vxm13 += 0.5*((dxmj/dx-0.5)^2) * (Vx[i,j  ] - 2.0*Vx[i+1,j  ] + Vx[i+2,j  ]);
                vxm24 += 0.5*((dxmj/dx-0.5)^2) * (Vx[i,j+1] - 2.0*Vx[i+1,j+1] + Vx[i+2,j+1]);
            end
        else
            if i>1
                vxm13 += 0.5*((dxmj/dx-0.5)^2) * (Vx[i-1,j  ] - 2.0*Vx[i,j  ] + Vx[i+1,j  ]);
                vxm24 += 0.5*((dxmj/dx-0.5)^2) * (Vx[i-1,j+1] - 2.0*Vx[i,j+1] + Vx[i+1,j+1]);
            end
        end
    end
    # Compute vx
    vxm = (1-dymi/dy) * vxm13 + (dymi/dy) * vxm24
    # println(p.x[k], " ", xv[k], " ", vxm)
    return vxm
end
export VxFromVxNodes

@views function VyFromVxNodes(Vy::Matrix{Float64}, k::Int64, p::Markers, xce::LinRange{Float64}, yv::LinRange{Float64}, dx::Float64, dy::Float64, ncx::Int64, ncy::Int64, new::Int64)
    # Interpolate vy
    i = Int64(round(trunc( (p.x[k] - xce[1])/dx ) + 1));
    j = Int64(round(trunc( (p.y[k] -  yv[1])/dy ) + 1));
    if i<1
        i=1
    elseif i>ncx+1
        i=ncx+1
    end
    if j<1
        j=1
    elseif j>ncy
        j = ncy
    end
    # Compute distances
    dxmj = p.x[k] - xce[i]
    dymi = p.y[k] -  yv[j]
    # Compute vy velocity for the left and right of the cell
    vym12 = Vy[i,j  ]*(1-dymi/dy) + Vy[i  ,j+1]*dymi/dy
    vym34 = Vy[i+1,j]*(1-dymi/dy) + Vy[i+1,j+1]*dymi/dy
    if new==1 
        if dymi/dy>=0.5
            if j<ncy
                vym12 += 0.5*((dymi/dy-0.5)^2) * ( Vy[i,j  ] - 2.0*Vy[i,j+1  ] + Vy[i,j+2  ]);
                vym34 += 0.5*((dymi/dy-0.5)^2) * ( Vy[i+1,j] - 2.0*Vy[i+1,j+1] + Vy[i+1,j+2]);
            end      
        else
            if j>1
                vym12 += 0.5*((dymi/dy-0.5)^2) * ( Vy[i,j-1  ] - 2.0*Vy[i,j  ] + Vy[i,j+1  ]);
                vym34 += 0.5*((dymi/dy-0.5)^2) * ( Vy[i+1,j-1] - 2.0*Vy[i+1,j] + Vy[i+1,j+1]);
            end
        end
    end
    # Compute vy
    vym = (1-dxmj/dx)*vym12 + (dxmj/dx)*vym34
    return vym
end
export VyFromVxNodes

@views function LocateMarkers(p::Markers,dx::Float64,dy::Float64,xc,yc,xmin::Float64,xmax::Float64,ymin::Float64,ymax::Float64)
    # Find marker cell indices
    @threads for k=1:p.nmark
        if (p.x[k]<xmin || p.x[k]>xmax || p.y[k]<ymin || p.y[k]>ymax) 
            p.phase[k] = -1
        end
        if p.phase[k]>=0
            dstx         = p.x[k] - xc[1]
            i            = ceil(Int, dstx / dx + 0.5)
            dsty         = p.y[k] - yc[1]
            j            = ceil(Int, dsty / dy + 0.5)
            p.cellx[k]   = i
            p.celly[k]   = j
        end
    end
end
export LocateMarkers

@views function Markers2Cells3!(p::Markers,phase::Matrix{Float64},xc::LinRange{Float64},yc::LinRange{Float64},dx::Float64,dy::Float64,ncx::Int64,ncy::Int64,prop::Vector{Float64},avg::Int64,mode::Int64)
    # if length(prop)  > p.nmark  println(length(prop)); println(p.nmark ); error("wtf") end
    phase0      = zeros(Float64, (ncx, ncy))
    phase0     .= phase
    weight      = zeros(Float64, (ncx, ncy))
    phase_th    = [similar(phase) for _ = 1:nthreads()] # per thread
    weight_th   = [similar(weight) for _ = 1:nthreads()] # per thread
    @threads for tid=1:nthreads()
        fill!(phase_th[tid] , 0)
        fill!(weight_th[tid], 0)
    end
    chunks = Iterators.partition(1:p.nmark, p.nmark ÷ nthreads())
    @sync for chunk in chunks
        Threads.@spawn begin
            tid = threadid()
            # fill!(phase_th[tid], 0)  # DON'T
            # fill!(weight_th[tid], 0)
            for k in chunk
                if p.phase[k]>=0
                # Get the indices:
                i = p.cellx[k]
                j = p.celly[k]
                # Relative distances
                dxm = 2.0 * abs(xc[i] - p.x[k])
                dym = 2.0 * abs(yc[j] - p.y[k])
                # Increment cell counts
                area = (1.0 - dxm / dx) * (1.0 - dym / dy)
                if mode==0 val  =  prop[Int64(p.phase[k])] end
                if mode==1 val  =  prop[k]                 end
                if avg==0 phase_th[tid][i,  j] += val       * area end
                if avg==1 phase_th[tid][i,  j] += (1.0/val) * area end
                if avg==2 phase_th[tid][i,  j] += log(val) * area end
                weight_th[tid][i, j] += area
                end
            end
        end
    end
    phase  .= reduce(+, phase_th)
    weight .= reduce(+, weight_th)
    phase ./= weight
    if avg==1
        phase .= 1.0 ./ phase
    end
    if avg==2
        phase .= exp.(phase)
    end
    phase[weight.<1e-10] .= phase0[weight.<1e-10]
    return
end
export Markers2Cells3!

@views function CountMarkers!(p::Markers,mpc,mpc_th,nmark,xmin::Float64,xmax::Float64,ymin::Float64,ymax::Float64,xc,yc,dx,dy::Float64,ncx::Int64,ncy::Int64)
    # Disable markers outside of the domain
    @threads for k=1:nmark
        if (p.x[k]<xmin || p.x[k]>xmax || p.y[k]<ymin || p.y[k]>ymax) 
            @inbounds p.phase[k] = -1
        end
    end

    # How many are outside? save indices for reuse
    nmark_out_th = zeros(Int64, nthreads())
    @threads for k=1:nmark
        if p.phase[k] == -1
            nmark_out_th[threadid()] += 1
        end
    end
    nmark_out = 0
    for ith=1:nthreads()
        nmark_out += nmark_out_th[ith]
    end
    @printf("%d markers out\n", nmark_out)

    # Count number of marker per cell
    @threads for j=1:ncy
        for i=1:ncx
            for ith=1:nthreads()
                mpc_th[ith,i,j] = 0.0
            end
        end
    end

    @threads for k=1:nmark # @avx ne marche pas ici
        if (p.phase[k]>=0)
            # Get the column:
            dstx = p.x[k] - xc[1];
            i    = Int64(round(ceil( (dstx/dx) + 0.5)));
            # Get the line:
            dsty = p.y[k] - yc[1];
            j    = Int64(round(ceil( (dsty/dy) + 0.5)));
            # Increment cell count
            mpc_th[threadid(),i,j] += 1.0
        end
    end

    @threads for j=1:ncy
        for i=1:ncx
            for ith=1:nthreads()
                if ith == 1 
                    mpc[i,j] = 0.0
                end
                mpc[i,j] += mpc_th[ith,i,j]
            end
        end
    end
    return       
end
export CountMarkers!

@views function CountMarkers3!( p::Markers, f::Fields2D, mat::MaterialParameters, dom::ModelDomain )
        
    ncx, ncy, dx, dy = dom.ncx, dom.ncy, dom.dx, dom.dy
    nphase = mat.nphase
    nmark  = p.nmark

    # Disable markers outside of the domain
    @threads for k=1:nmark
        if (p.x[k]<dom.xmin || p.x[k]>dom.xmax || p.y[k]<dom.ymin || p.y[k]>dom.ymax) 
            @inbounds p.phase[k] = -1
        end
    end

    # How many are outside? save indices for reuse
    nmark_out_th = zeros(Int64, nthreads())
    @threads for k=1:nmark
        if p.phase[k] == -1
            nmark_out_th[threadid()] += 1
        end
    end
    nmark_out = 0
    for ith=1:nthreads()
        nmark_out += nmark_out_th[ith]
    end
    @printf("%d markers out\n", nmark_out)

    # Count number of marker per cell
    @threads for j=1:ncy
        for i=1:ncx
            for ith=1:nthreads()
                f.mpc_th[ith,i,j] = 0.0
                for m=1:nphase
                    f.phase_perc_th[ith,m,i,j] = 0.0
                end
            end
        end
    end

    @threads for k=1:nmark # @avx ne marche pas ici
        if (p.phase[k]>=0)
            m = Int64(p.phase[k])
            # Get the column:
            dstx = p.x[k] - dom.xc[1];
            i    = Int64(round(ceil( (dstx/dx) + 0.5)));
            # Get the line:
            dsty = p.y[k] - dom.yc[1];
            j    = Int64(round(ceil( (dsty/dy) + 0.5)));
            # Increment cell count
            f.mpc_th[threadid(),i,j] += 1.0
            f.phase_perc_th[threadid(),m,i,j] += 1.0
        end
    end

    @threads for j=1:ncy
        for i=1:ncx
            for ith=1:nthreads()
                if ith == 1 
                    f.mpc[i,j] = 0.0
                    for m=1:nphase
                        f.phase_perc[m,i,j] = 0.0
                    end
                end
                f.mpc[i,j] += f.mpc_th[ith,i,j]
                for m=1:nphase
                    f.phase_perc[m,i,j] += f.phase_perc_th[ith,m,i,j]
                end
            end
        end
    end

    # Normalize phase percentages
    for m=1:nphase
        f.phase_perc[m,:,:] ./= f.mpc
    end
      
    # Show infos
    println("markers per cell:")
    println("min: ", minimum(f.mpc), " --- max: ", maximum(f.mpc))
    for m=1:nphase
        println("phase ", m, " :")
        println("min: ", minimum(f.phase_perc[m,:,:]), " --- max: ", maximum(f.phase_perc[m,:,:]))
    end
end
export CountMarkers3!

@views function ReseedMarkers2!( p::Markers, mpc::Matrix{Float64}, dom::ModelDomain)

    ncx, ncy, dx, dy = dom.ncx, dom.ncy, dom.dx, dom.dy
    mlist       = Array{Array{Int64}}(undef,ncx,ncy)
    mlist_th    = Array{Array{Int64}}(undef,nthreads(),ncx,ncy)

    for i=1:ncx
        for j=1:ncy
            mlist[i,j] = zeros(Int64(mpc[i,j]))
        end
    end

    count = zeros(Int64,ncx,ncy)

    for k=1:p.nmark # @avx ne marche pas ici
        if p.phase[k]>=0
            # Get the column:
            dstx = p.x[k] - dom.xc[1];
            i    = Int64(round(ceil( (dstx/dx) + 0.5)));
            # Get the line:
            dsty = p.y[k] - dom.yc[1];
            j    = Int64(round(ceil( (dsty/dy) + 0.5)));
            # Increment cell count
            count[i,j] += 1
            mlist[i,j][count[i,j]] = k
        end
    end

    # display(mlist)
    dxy = 10*(dx*dx + dy*dy)
    for i=1:ncx
        for j=1:ncy

            if mpc[i,j] < 4

                new_x = dom.xc[i]-dx/4.0
                new_y = dom.yc[j]-dy/4.0
                ineigh = FindClosestNeighbour( p, mpc, mlist, new_x, new_y, dxy, i, j, ncx, ncy )
                if ineigh>0
                    p.nmark+=1
                    p.x[p.nmark] = new_x
                    p.y[p.nmark] = new_y
                    CopyMarkers!(p,ineigh)
                end

                new_x = dom.xc[i]+dx/4.0
                new_y = dom.yc[j]-dy/4.0
                ineigh = FindClosestNeighbour( p, mpc, mlist, new_x, new_y, dxy, i, j, ncx, ncy )
                if ineigh>0
                    p.nmark+=1
                    p.x[p.nmark] = new_x
                    p.y[p.nmark] = new_y
                    CopyMarkers!(p,ineigh)
                end

                new_x = dom.xc[i]+dx/4.0
                new_y = dom.yc[j]+dy/4.0
                ineigh = FindClosestNeighbour( p, mpc, mlist, new_x, new_y, dxy, i, j, ncx, ncy )
                if ineigh>0
                    p.nmark+=1
                    p.x[p.nmark] = new_x
                    p.y[p.nmark] = new_y
                    CopyMarkers!(p,ineigh)
                end

                new_x = dom.xc[i]-dx/4.0
                new_y = dom.yc[j]+dy/4.0
                ineigh = FindClosestNeighbour( p, mpc, mlist, new_x, new_y, dxy, i, j, ncx, ncy )
                if ineigh>0
                    p.nmark+=1
                    p.x[p.nmark] = new_x
                    p.y[p.nmark] = new_y
                    CopyMarkers!(p,ineigh)
                end

            end
        end
    end
end
export ReseedMarkers2!

@views function LocateMarkers2!(p::Markers,dom::ModelDomain)
    # Find marker cell indices
    @threads for k=1:p.nmark
        if (p.x[k]<dom.xmin || p.x[k]>dom.xmax || p.y[k]<dom.ymin || p.y[k]>dom.ymax) 
            p.phase[k] = -1
        end
        if p.phase[k]>=0
            dstx         = p.x[k] - dom.xc[1]
            i            = ceil(Int, dstx / dom.dx + 0.5)
            dsty         = p.y[k] - dom.yc[1]
            j            = ceil(Int, dsty / dom.dy + 0.5)
            p.cellx[k]   = i
            p.celly[k]   = j
        end
    end
end
export LocateMarkers2!

################################################################
################################################################

function RungeKutta2!(p::Markers, f::Fields2D, params::ModelParameters, BC::BoundaryConditions, dom::ModelDomain)
    
    nmark = p.nmark
    rkv, rkw, dt   = params.rkv, params.rkw, params.dt
    ncx, ncy, dx, dy = dom.ncx, dom.ncy, dom.dx, dom.dy
 
    # Marker advection with 4th order Roger-Kutta
    @threads for k=1:nmark
        in = p.phase[k]>=0
            x0 = p.x[k];
            y0 = p.y[k];
            vx = 0.0
            vy = 0.0
            # Roger-Gunther loop
            for rk=1:4
                # println("step", rk)
                # Interp velocity from grid
                vxm = VxFromVxNodes(f.Vx, k, p, dom.xv, dom.yce, dx, dy, ncx, ncy, 1)
                vym = VyFromVxNodes(f.Vy, k, p, dom.xce, dom.yv, dx, dy, ncx, ncy, 1)
                # Temporary RK advection steps
                p.x[k] = x0 + rkv[rk]*dt*vxm
                p.y[k] = y0 + rkv[rk]*dt*vym
                if BC.periodix==1
                    if p.x[k]<dom.xv[1]
                        # println("periW 1 ", p.x[k], vxm)
                        p.x[k]+=(dom.xv[end]-dom.xv[1]) 
                        # println("periW 2 ", p.x[k], vxm)
                    end
                    if p.x[k]>dom.xv[end] 
                        # println("periE")
                        p.x[k]-=(dom.xv[end]-dom.xv[1]) 
                    end
                end
                # Average final velocity 
                vx    += rkw[rk]*vxm
                vy    += rkw[rk]*vym
            end
            # Advect points
            p.x[k] = x0 + (in==1) * rkv[4]*dt*vx
            p.y[k] = y0 + (in==1) * rkv[4]*dt*vy
            if BC.periodix==1
                if p.x[k]<dom.xv[1]
                    # println("periW 1 ", p.x[k], vxm)
                    p.x[k]+=(dom.xv[end]-dom.xv[1]) 
                    # println("periW 2 ", p.x[k], vxm)
                end
                if p.x[k]>dom.xv[end] 
                    # println("periE")
                    p.x[k]-=(dom.xv[end]-dom.xv[1]) 
                end
            end
    end
    end
    export RungeKutta2!
