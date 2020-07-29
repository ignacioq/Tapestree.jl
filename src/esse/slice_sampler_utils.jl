#=

Slice sampling utilities

Ignacio Quintero Mächler

t(-_-t)

September 23 2017

=#



"""
    w_sampler(lhf         ::Function, 
              p           ::Array{Array{Float64,1},1},
              fp          ::Array{Array{Float64,1},1},
              nnps        ::Array{Int64,1},
              nps         ::Array{Int64,1},
              phid        ::Array{Int64,1},
              mvps        ::Array{Array{Int64,1},1},
              nngps        ::Array{Array{Bool,1},1},
              mvhfs       ::Array{Array{Int64,1},1},
              hfgps       ::Array{Array{Bool,1},1},
              npars       ::Int64,
              optimal_w   ::Float64,
              screen_print::Int64,
              nburn       ::Int64,
              ntakew      ::Int64,
              nswap       ::Int64,
              ncch        ::Int64,
              winit       ::Float64)

Run slice sampler for burn-in and to estimate appropriate w's.
"""
function w_sampler(lhf         ::Function, 
                   p           ::DArray{Array{Float64,1},1,Array{Array{Float64,1},1}},
                   fp          ::DArray{Array{Float64,1},1,Array{Array{Float64,1},1}},
                   nnps        ::Array{Int64,1},
                   nps         ::Array{Int64,1},
                   phid        ::Array{Int64,1},
                   mvps        ::Array{Array{Int64,1},1},
                   nngps        ::Array{Array{Bool,1},1},
                   mvhfs       ::Array{Array{Int64,1},1},
                   hfgps       ::Array{Array{Bool,1},1},
                   npars       ::Int64,
                   optimal_w   ::Float64,
                   screen_print::Int64,
                   nburn       ::Int64,
                   ntakew      ::Int64,
                   nswap       ::Int64,
                   ncch        ::Int64,
                   winit       ::Float64,
                   dt          ::Float64)

  if nburn < ntakew
    ntakew = nburn
  end

  # maximum number of multivariate updates
  maxmvu = if iszero(lastindex(mvps)) && iszero(lastindex(mvhfs))
    zero(1)
  elseif iszero(lastindex(mvhfs))
    maximum(map(length,mvps))
  elseif iszero(lastindex(mvps))
    maximum(map(length,mvhfs))
  else
    maximum((maximum(map(length,mvps)), maximum(map(length,mvhfs))))
  end

  # temperature
  t, o = make_temperature(dt, ncch)

  Lv = Array{Float64,1}(undef, maxmvu)
  Rv = Array{Float64,1}(undef, maxmvu)

  # preallocate pp and fpp
  pp  = copy(p[1])
  fpp = copy(fp[1])

  w  = fill(winit, npars)

  # Make distributed array for parameters
  ps = [Array{Float64,2}(undef, nburn, npars) for i in Base.OneTo(ncch)]
  ps = distribute(ps)

  # Make distributed array for iterations per chain `ipc`
  ipc = SharedArray{Int64,1}(zeros(Int64,ncch))

  # make SharedArray for current posteriors
  lhc  = [lhf(p[c], fp[c], t[o[c]]) for c in Base.OneTo(ncch)]
  lhc = SharedArray(lhc)

  # make SharedArray for temperature order `o`
  o = SharedArray(o)

  # pre-estimate chain swaps
  tns   = fld(nburn,nswap)
  allsw = [wchains(o) for i in Base.OneTo(tns)]
  ssper = SharedArray{Bool,1}(fill(false, tns))
  sspee = SharedArray{Bool,1}(fill(false, tns))

  # length fp
  lfp = length(fp[1])


  @sync @distributed for c in Base.OneTo(ncch)

    # swap variables
    lswap, ws = 0, 0

    for it in Base.OneTo(nburn)

      # slice parameter cycle
      lhc[c] = 
        slice_cycle(lhf, lhc[c], p[c], fp[c], pp, fpp, Lv, Rv, 
          nnps, nps, phid, mvps, nngps, mvhfs, hfgps, w, npars, lfp, t[o[c]])

      #log chain current iter
      ipc[c] = it

      # log parameters
      @inbounds ps[c][it, :] = p[c]

      ## swap chains
      lswap += 1
      if lswap == nswap
        ws += 1

        j, k = allsw[ws]

        cij = c == j   # c is j
        cik = c == k   # c is k

        println("swap for j = ", j, " <-> k = ",k)
        # swapper
        if cij
          println("swapper came in swap")
          # wait for other chain
          while true
            ipc[c] == ipc[k] && break
            #println("iter c:", ipc[c], "; iter k:", ipc[k])
          end
          println(c, " was able to move onto swap ", ipc)
          # make swap
          ssper[ws] = swap_chains!(j, k, o, t, lhc)
          println("SWAP for ", j," ", k)
          while true 
            sspee[ws] && break
          end
        end

        # swappee
        if cik
          println("swappee came in swap")
          # wait for other chain and swap to complete
          while true
            ipc[c] == ipc[j] && ssper[ws] && break
            #println("iter c:", ipc[c], "; iter j:", ipc[j])
          end
          sspee[ws] = true
          println("swappee moved on from swap")
        end

        lswap = 0
        println(c, " moved out of swap")
      end
    end
  end



  sps = nburn-ntakew

  ps = ps[findfirst(x -> isone(x), o)][(nburn-ntakew+1):nburn,:]

  w = optimal_w .* (reduce(max, ps, dims=1) .- reduce(min, ps, dims=1))
  w = reshape(w, size(w,2))

  return p, fp, w, o, t
end







"""
    loop_slice_sampler(lhf         ::Function, 
                       p           ::Array{Array{Float64,1},1},
                       fp          ::Array{Array{Float64,1},1},
                       nnps        ::Array{Int64,1},
                       nps         ::Array{Int64,1},
                       phid        ::Array{Int64,1},
                       mvps        ::Array{Array{Int64,1},1},
                       nngps       ::Array{Array{Bool,1},1},
                       mvhfs       ::Array{Array{Int64,1},1},
                       hfgps       ::Array{Array{Bool,1},1},
                       w           ::Array{Float64,1},
                       npars       ::Int64,
                       niter       ::Int64,
                       nthin       ::Int64,
                       nswap       ::Int64,
                       ncch        ::Int64,
                       o          ::Array{Int64,1}, 
                       t          ::Array{Float64,1},
                       screen_print::Int64)

Run slice-sampling.
"""
function loop_slice_sampler(lhf         ::Function, 
                            p           ::Array{Array{Float64,1},1},
                            fp          ::Array{Array{Float64,1},1},
                            nnps        ::Array{Int64,1},
                            nps         ::Array{Int64,1},
                            phid        ::Array{Int64,1},
                            mvps        ::Array{Array{Int64,1},1},
                            nngps       ::Array{Array{Bool,1},1},
                            mvhfs       ::Array{Array{Int64,1},1},
                            hfgps       ::Array{Array{Bool,1},1},
                            w           ::Array{Float64,1},
                            npars       ::Int64,
                            niter       ::Int64,
                            nthin       ::Int64,
                            nswap       ::Int64,
                            ncch        ::Int64,
                            o          ::Array{Int64,1}, 
                            t          ::Array{Float64,1},
                            screen_print::Int64)

  # maximum number of parameters in multivariate updates
  maxmvu = if iszero(lastindex(mvps)) && iszero(lastindex(mvhfs))
    zero(1)
  elseif iszero(lastindex(mvhfs))
    maximum(map(length,mvps))
  elseif iszero(lastindex(mvps))
    maximum(map(length,mvhfs))
  else
    maximum((maximum(map(length,mvps)), maximum(map(length,mvhfs))))
  end

  nlogs = fld(niter,nthin)

  #preallocate logging arrays
  its  =  Array{Float64,1}(undef, nlogs)
  hlog =  Array{Float64,2}(undef, nlogs, ncch)
  ps   = [Array{Float64,2}(undef, nlogs, npars) for i in Base.OneTo(ncch)]

  # preallocate chains order
  olog   = Array{Int64,2}(undef, nlogs, ncch)

  # preallocate changing vectors
  Lv    = Array{Float64,1}(undef, maxmvu)
  Rv    = Array{Float64,1}(undef, maxmvu)

  lthin, lit, lswap = 0, 0, 0

  # preallocate pp and fpp
  pp  = copy(p[1])
  fpp = copy(fp[1])

  # length fp
  lfp = length(fp[1])

  # start iterations
  prog = Progress(niter, screen_print, "running slice-sampler...", 20)

  # starting posteriors
  lhc = [lhf(p[c], fp[c], t[o[c]]) for c in Base.OneTo(ncch)]

  for it in Base.OneTo(niter) 
    for c in Base.OneTo(ncch)
      lhc[c] = 
        slice_cycle(lhf, lhc[c], p[c], fp[c], pp, fpp,
          nnps, nps, phid, mvps, nngps, mvhfs, hfgps, w, npars, t[o[c]])
    end

    # log samples
    lthin += 1
    if lthin == nthin
      @inbounds begin
        lit += 1
        setindex!(its,  it,  lit)
        setindex!(hlog, lhc, lit, :)
        for c in Base.OneTo(ncch)
          setindex!(ps[c], p[c], lit, :)
        end
        setindex!(olog, o, lit, :)
      end
      lthin = 0
    end

    # swap chains
    if ncch > 1
      lswap += 1
      if lswap == nswap
        swap_chains!(o, t, lhc)
        lswap = 0
      end
    end

    next!(prog)
  end

  return its, hlog, ps, olog
end





"""
    slice_cycle(lhf  ::Function, 
                lhc  ::Float64,
                p    ::Array{Float64,1},
                fp   ::Array{Float64,1},
                pp   ::Array{Float64,1}, 
                fpp  ::Array{Float64,1}, 
                Lv   ::Array{Float64,1},
                Rv   ::Array{Float64,1},
                nnps ::Array{Int64,1},
                nps  ::Array{Int64,1},
                phid ::Array{Int64,1},
                mvps ::Array{Array{Int64,1},1},
                nngps::Array{Array{Bool,1},1},
                mvhfs::Array{Array{Int64,1},1},
                hfgps::Array{Array{Bool,1},1},
                w    ::Array{Float64,1},
                npars::Int64,
                lfp  ::Int64,
                ti   ::Float64)

A full parameter cycle for slice-sampling.
"""
function slice_cycle(lhf  ::Function, 
                     lhc  ::Float64,
                     p    ::Array{Float64,1},
                     fp   ::Array{Float64,1},
                     pp   ::Array{Float64,1}, 
                     fpp  ::Array{Float64,1}, 
                     Lv   ::Array{Float64,1},
                     Rv   ::Array{Float64,1},
                     nnps ::Array{Int64,1},
                     nps  ::Array{Int64,1},
                     phid ::Array{Int64,1},
                     mvps ::Array{Array{Int64,1},1},
                     nngps::Array{Array{Bool,1},1},
                     mvhfs::Array{Array{Int64,1},1},
                     hfgps::Array{Array{Bool,1},1},
                     w    ::Array{Float64,1},
                     npars::Int64,
                     lfp  ::Int64,
                     ti   ::Float64)

  #=
  univariate updates
  =#
  # nonnegative parameters
  for j in nnps
    S    = lhc - randexp()
    L, R = find_nonneg_int(p, pp, fp, j, S, lhf, w[j], ti, npars)
    lhc  = sample_int(p, pp, fp, j, L, R, S, lhf, ti, npars)
  end

  # real line parameters
  for j in nps
    S    = lhc - randexp()
    L, R = find_real_int(p, pp, fp, j, S, lhf, w[j], ti, npars)
    lhc  = sample_int(p, pp, fp, j, L, R, S, lhf, ti, npars)
  end

  # hidden factors
  for j in phid
    S    = lhc - randexp()
    L, R = 
      find_nonneg_int(p, pp, fp, fpp, j, S, lhf, w[j], ti, npars, lfp)
    lhc = 
      sample_int(p, pp, fp, fpp, j, L, R, S, lhf, ti, npars, lfp)
  end

  #=
  multivariate updates
  =#
  # non hidden factors
  for (i, mvp) in enumerate(mvps)
    S = lhc - randexp()
    @views find_rect(p, pp, fp, mvp, nngps[i], S, lhf, w, Lv, Rv, ti, npars)
    lhc = sample_rect(p, pp, fp, Lv, Rv, mvp, S, lhf, ti,  npars)
  end

  # for joint hidden factors
  for (i, mvp) in enumerate(mvhfs)
    S = lhc - randexp()
    @views find_rect(p, pp, fp, fpp, mvp, hfgps[i], S, lhf, w, Lv, Rv, ti, npars, lfp)
    @views lhc = 
      sample_rect(p, pp, fp, fpp, Lv, Rv, mvp, hfgps[i], S, lhf, ti, npars, lfp)
  end

  return lhc
end








"""
    find_nonneg_int(p    ::Array{Float64}, 
                    pp   ::Array{Float64},
                    j    ::Int64, 
                    S    ::Float64, 
                    postf::Function, 
                    w    ::Float64)

Estimate a nonnegative slice interval.
"""
function find_nonneg_int(p    ::Array{Float64,1}, 
                         pp   ::Array{Float64,1},
                         fp   ::Array{Float64,1},
                         j    ::Int64, 
                         S    ::Float64, 
                         lhf::Function, 
                         w    ::Float64,
                         ti   ::Float64,
                         npars::Int64)

  unsafe_copyto!(pp, 1, p, 1, npars)

  L = pp[j] - w*rand()
  R = L + w

  if L <= 0.0
    L = 1e-30
  end

  # left extreme
  pp[j] = L
  while S < lhf(pp, fp, ti)
    L -= w
    if L <= 0.0
      L = 1e-30
      break
    end
    pp[j] = L
  end

  # right extreme
  pp[j] = R
  while S < lhf(pp, fp, ti)
    R    += w
    pp[j] = R
  end

  return (L, R)::NTuple{2,Float64}
end





"""
    find_nonneg_int(p    ::Array{Float64,1}, 
                    pp   ::Array{Float64,1},
                    fp   ::Array{Float64,1},
                    fpp  ::Array{Float64,1},
                    j    ::Int64, 
                    S    ::Float64, 
                    postf::Function, 
                    w    ::Float64)

Estimate a non_negative slice interval for hidden factors.
"""
function find_nonneg_int(p    ::Array{Float64,1}, 
                         pp   ::Array{Float64,1},
                         fp   ::Array{Float64,1},
                         fpp  ::Array{Float64,1},
                         j    ::Int64, 
                         S    ::Float64, 
                         lhf  ::Function, 
                         w    ::Float64,
                         ti    ::Float64,
                         npars::Int64,
                         lfp  ::Int64)

  unsafe_copyto!(pp,  1, p,  1, npars)
  unsafe_copyto!(fpp, 1, fp, 1, lfp)

  L = fpp[j] - w*rand()
  R = L + w

  if L <= 0.0
    L = 1e-30
  end

  # left extreme
  fpp[j] = L
  while S < lhf(pp, fpp, ti)
    L -= w
    if L <= 0.0
      L = 1e-30
      break
    end
    fpp[j] = L
  end

  # right extreme
  fpp[j] = R
  while S < lhf(pp, fpp, ti)
    R    += w
    fpp[j] = R
  end

  return (L, R)::NTuple{2,Float64}
end





"""
    find_real_int(p    ::Array{Float64}, 
                  pp   ::Array{Float64}, 
                  fp   ::Array{Float64},
                  j    ::Int64, 
                  S    ::Float64, 
                  postf::Function, 
                  w    ::Float64)

Estimate a nonnegative slice interval.
"""
function find_real_int(p    ::Array{Float64,1}, 
                       pp   ::Array{Float64,1}, 
                       fp   ::Array{Float64,1}, 
                       j    ::Int64, 
                       S    ::Float64, 
                       lhf  ::Function, 
                       w    ::Float64,
                       ti   ::Float64,
                       npars::Int64)

  unsafe_copyto!(pp, 1, p, 1, npars)

  L = pp[j] - w*rand()
  R = L + w

  # left extreme
  pp[j] = L::Float64
  while S < lhf(pp, fp, ti)
    L    -= w::Float64
    pp[j] = L::Float64
  end

  # right extreme
  pp[j] = R::Float64
  while S < lhf(pp, fp, ti)
    R    += w::Float64
    pp[j] = R::Float64
  end

  return (L, R)::NTuple{2,Float64}
end




"""
    sample_int(p    ::Array{Float64,1}, 
               pp   ::Array{Float64,1},
               fp   ::Array{Float64,1},
               j    ::Int64, 
               L    ::Float64, 
               R    ::Float64, 
               S    ::Float64, 
               postf::Function)

Take one sample within the interval of the slice.
"""
function sample_int(p    ::Array{Float64,1}, 
                    pp   ::Array{Float64,1},
                    fp   ::Array{Float64,1},
                    j    ::Int64, 
                    L    ::Float64, 
                    R    ::Float64, 
                    S    ::Float64, 
                    lhf  ::Function,
                    ti   ::Float64,
                    npars::Int64)

  @inbounds begin
    unsafe_copyto!(pp, 1, p, 1, npars)

    while true
      pp[j] = (L + rand()*(R-L))::Float64

      hc = lhf(pp, fp, ti)::Float64
      if S < hc
        unsafe_copyto!(p, 1, pp, 1, npars)
        return hc::Float64
      end

      if pp[j] < p[j]
        L = pp[j]::Float64
      else
        R = pp[j]::Float64
      end
    end

  end
end





"""
    sample_int(p    ::Array{Float64,1}, 
               pp   ::Array{Float64,1},
               fp   ::Array{Float64,1},
               fpp  ::Array{Float64,1},
               j    ::Int64, 
               L    ::Float64, 
               R    ::Float64, 
               S    ::Float64, 
               postf::Function)

Take one sample within the interval of the slice for hidden factors.
"""
function sample_int(p    ::Array{Float64,1}, 
                    pp   ::Array{Float64,1},
                    fp   ::Array{Float64,1},
                    fpp  ::Array{Float64,1},
                    j    ::Int64, 
                    L    ::Float64, 
                    R    ::Float64, 
                    S    ::Float64, 
                    lhf  ::Function, 
                    ti   ::Float64,
                    npars::Int64,
                    lfp  ::Int64)

  @inbounds begin
    unsafe_copyto!(pp,  1, p,  1, npars)
    unsafe_copyto!(fpp, 1, fp, 1, lfp)

    while true
      fpp[j] = (L + rand()*(R-L))::Float64

      hc = lhf(pp, fpp, ti)::Float64
      if S < hc
        unsafe_copyto!(p,  1, pp,  1, npars)
        unsafe_copyto!(fp, 1, fpp, 1, lfp)
        return hc::Float64
      end

      if fpp[j] < fp[j]
        L = fpp[j]::Float64
      else
        R = fpp[j]::Float64
      end
    end

  end
end





"""
    find_rect(p    ::Array{Float64,1}, 
              pp   ::Array{Float64,1}, 
              fp   ::Array{Float64,1}, 
              mvp  ::Array{Int64,1}, 
              nngp  ::Array{Bool,1},
              S    ::Float64, 
              postf::Function, 
              w    ::Array{Float64,1},
              Lv   ::Array{Float64,1},
              Rv   ::Array{Float64,1},
              npars::Int64)

Returns left and right vertices of an hyper-rectangle.
"""
function find_rect(p    ::Array{Float64,1}, 
                   pp   ::Array{Float64,1}, 
                   fp   ::Array{Float64,1}, 
                   mvp  ::Array{Int64,1}, 
                   nngp  ::Array{Bool,1},
                   S    ::Float64, 
                   lhf  ::Function, 
                   w    ::Array{Float64,1},
                   Lv   ::Array{Float64,1},
                   Rv   ::Array{Float64,1},
                   ti   ::Float64,
                   npars::Int64)

  @inbounds begin

    unsafe_copyto!(pp, 1, p, 1, npars)

    # randomly start rectangle
    for (i,j) in enumerate(mvp)
      Lv[i] = pp[j] - w[j]*rand()
      if nngp[i] && Lv[i] <= 0.0
        Lv[i] = 1e-30
      end
      Rv[i] = Lv[i] + w[j]
    end

    # left extremes
    for (i,j) in enumerate(mvp)
      pp[j] = Lv[i]::Float64
      while S < lhf(pp, fp, ti)
        Lv[i] -= w[j]::Float64
        if nngp[i] && Lv[i] <= 0.0
          Lv[i] = 1e-30
          break
        end
        pp[j]  = Lv[i]::Float64
      end
    end

    # right extremes
    for (i,j) in enumerate(mvp)
      pp[j] = Rv[i]::Float64
      while S < lhf(pp, fp, ti)
        Rv[i] += w[j]::Float64
        pp[j]  = Rv[i]::Float64
      end
    end
  end

  return nothing
end





"""
    find_rect(p    ::Array{Float64,1}, 
              pp   ::Array{Float64,1}, 
              fp   ::Array{Float64,1}, 
              fpp  ::Array{Float64,1},
              mvp  ::Array{Int64,1}, 
              hfgp ::Array{Bool,1},
              S    ::Float64, 
              postf::Function, 
              w    ::Array{Float64,1},
              Lv   ::Array{Float64,1},
              Rv   ::Array{Float64,1},
              npars::Int64)

Returns left and right vertices of an hyper-rectangle for mixed hidden factors.
"""
function find_rect(p    ::Array{Float64,1}, 
                   pp   ::Array{Float64,1}, 
                   fp   ::Array{Float64,1}, 
                   fpp  ::Array{Float64,1},
                   mvp  ::Array{Int64,1}, 
                   hfgp ::Array{Bool,1},
                   S    ::Float64, 
                   lhf  ::Function, 
                   w    ::Array{Float64,1},
                   Lv   ::Array{Float64,1},
                   Rv   ::Array{Float64,1},
                   ti   ::Float64,
                   npars::Int64,
                   lfp  ::Int64)

  @inbounds begin

    unsafe_copyto!(pp,  1, p,  1, npars)
    unsafe_copyto!(fpp, 1, fp, 1, lfp)

    # randomly start rectangle
    for (i,j) in enumerate(mvp)
      if hfgp[i]
        Lv[i] = fpp[j] - w[j]*rand()
        if Lv[i] <= 0.0
          Lv[i] = 1e-30
        end
        Rv[i] = Lv[i] + w[j]
      else
        Lv[i] = pp[j] - w[j]*rand()
        Rv[i] = Lv[i] + w[j]
      end
    end

    # left extremes
    for (i,j) in enumerate(mvp)
      if hfgp[i]
        fpp[j] = Lv[i]::Float64
        while S < lhf(pp, fpp, ti)
          Lv[i] -= w[j]::Float64
          if Lv[i] <= 0.0
            Lv[i] = 1e-30
            break
          end
          fpp[j] = Lv[i]::Float64
        end
      else
        pp[j] = Lv[i]::Float64
        while S < lhf(pp, fpp, ti)
          Lv[i] -= w[j]::Float64
          pp[j]  = Lv[i]::Float64
        end
      end
    end

    # right extremes
    for (i,j) in enumerate(mvp)
      if hfgp[i]
        fpp[j] = Rv[i]::Float64
        while S < lhf(pp, fpp, ti)
          Rv[i] += w[j]::Float64
          fpp[j]  = Rv[i]::Float64
        end
      else
        pp[j] = Rv[i]::Float64
        while S < lhf(pp, fpp, ti)
          Rv[i] += w[j]::Float64
          pp[j]  = Rv[i]::Float64
        end
      end
    end
  end

  return nothing
end





"""
    sample_rect(p    ::Array{Float64,1}, 
                pp   ::Array{Float64,1},
                fp   ::Array{Float64,1},
                Lv   ::Array{Float64,1}, 
                Rv   ::Array{Float64,1}, 
                mvp  ::Array{Float64,1},
                S    ::Float64, 
                postf::Function, 
                npars::Int64)

Sample from an hyper-rectangle with shrinking.
"""
function sample_rect(p    ::Array{Float64,1}, 
                     pp   ::Array{Float64,1},
                     fp   ::Array{Float64,1},
                     Lv   ::Array{Float64,1}, 
                     Rv   ::Array{Float64,1}, 
                     mvp  ::Array{Int64,1},
                     S    ::Float64, 
                     lhf  ::Function,
                     ti   ::Float64,
                     npars::Int64)

  @inbounds begin

    unsafe_copyto!(pp, 1, p, 1, npars)

    while true

      # multivariate uniform sample
      for (i,j) in enumerate(mvp)
        pp[j] = Lv[i] + rand()*(Rv[i] - Lv[i])
      end

      # posterior
      hc = lhf(pp, fp, ti)

      if S < hc
        unsafe_copyto!(p, 1, pp, 1, npars)
        return hc::Float64
      end

      for (i,j) in enumerate(mvp)
        if pp[j] < p[j]
            Lv[i] = pp[j]::Float64
        else
            Rv[i] = pp[j]::Float64
        end
      end

    end
  end

  return nothing
end





"""
    sample_rect(p    ::Array{Float64,1}, 
                pp   ::Array{Float64,1},
                fp   ::Array{Float64,1},
                fpp  ::Array{Float64,1},
                Lv   ::Array{Float64,1}, 
                Rv   ::Array{Float64,1}, 
                mvp  ::Array{Int64,1},
                S    ::Float64, 
                postf::Function, 
                npars::Int64)

Sample from an hyper-rectangle with shrinking for hidden factors.
"""
function sample_rect(p    ::Array{Float64,1}, 
                     pp   ::Array{Float64,1},
                     fp   ::Array{Float64,1},
                     fpp  ::Array{Float64,1},
                     Lv   ::Array{Float64,1}, 
                     Rv   ::Array{Float64,1}, 
                     mvp  ::Array{Int64,1},
                     hfgp ::Array{Bool,1},
                     S    ::Float64, 
                     lhf  ::Function, 
                     ti   ::Float64,
                     npars::Int64,
                     lfp  ::Int64)

  @inbounds begin

    unsafe_copyto!(pp,  1, p,  1, npars)
    unsafe_copyto!(fpp, 1, fp, 1, lfp)

    while true

      # multivariate uniform sample
      for (i,j) in enumerate(mvp)
        if hfgp[i]
          fpp[j] = Lv[i] + rand()*(Rv[i] - Lv[i])
        else
          pp[j] = Lv[i] + rand()*(Rv[i] - Lv[i])
        end
      end

      # posterior
      hc = lhf(pp, fpp, ti)

      if S < hc
        unsafe_copyto!(p, 1, pp, 1, npars)
        unsafe_copyto!(fp, 1, fpp, 1, lfp)
        return hc::Float64
      end

      for (i,j) in enumerate(mvp)
        if hfgp[i]
          if fpp[j] < fp[j]
              Lv[i] = fpp[j]::Float64
          else
              Rv[i] = fpp[j]::Float64
          end
        else
         if pp[j] < p[j]
              Lv[i] = pp[j]::Float64
          else
              Rv[i] = pp[j]::Float64
          end
        end
      end

    end
  end

  return nothing
end






"""
    scaleT(T::Float64, rate::Float64)

Make scaling function given the objective acceptance rates.
"""
function scaleT(T::Float64, rate::Float64)
  if 0.4 > rate
    T /= (2.0 - rate/0.4)
  else
    T *= (1.0 + (rate - 0.4)/0.6)
  end

  return T
end



