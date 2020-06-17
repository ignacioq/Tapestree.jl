#=

Slice sampling utilities

Ignacio Quintero Mächler

t(-_-t)

September 23 2017

=#





"""
    loop_slice_sampler(lhf         ::Function, 
                       p           ::Array{Float64,1},
                       fp          ::Array{Float64,1},
                       nnps        ::Array{Int64,1},
                       nps         ::Array{Int64,1},
                       phid        ::Array{Int64,1},
                       mvps        ::Array{Array{Int64,1},1},
                       mvhfs       ::Array{Array{Int64,1},1},
                       w           ::Array{Float64,1},
                       npars       ::Int64,
                       niter       ::Int64,
                       nthin       ::Int64,
                       screen_print::Int64)

Run slice sampling.
"""
function loop_slice_sampler(lhf         ::Function, 
                            p           ::Array{Float64,1},
                            fp          ::Array{Float64,1},
                            nnps        ::Array{Int64,1},
                            nps         ::Array{Int64,1},
                            phid        ::Array{Int64,1},
                            mvps        ::Array{Array{Int64,1},1},
                            nngps       ::Array{Array{Bool,1},1},
                            mvhfs       ::Array{Array{Int64,1},1},
                            w           ::Array{Float64,1},
                            npars       ::Int64,
                            niter       ::Int64,
                            nthin       ::Int64,
                            screen_print::Int64)

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

  nlogs = fld(niter,nthin)
  its   = Array{Float64,1}(undef, nlogs)
  hlog  = Array{Float64,1}(undef, nlogs)
  ps    = Array{Float64,2}(undef, nlogs,npars)
  Lv    = Array{Float64,1}(undef, maxmvu)
  Rv    = Array{Float64,1}(undef, maxmvu)

  lthin, lit = 0, 0

  # preallocate pp
  pp = copy(p)

  # preallocate fpp
  fpp = copy(fp)

  # start iterations
  prog = Progress(niter, screen_print, "running slice-sampler...", 20)

  hc = lhf(p, fp)

  for it in Base.OneTo(niter) 

    #=
    univariate updates
    =#
    # nonnegative parameters
    for j in nnps
     S     = hc - randexp()
     L, R  = find_nonneg_int(p, pp, fp, j, S, lhf, w[j], npars)
     p, hc = sample_int(p, pp, fp, j, L, R, S, lhf, npars)
    end

    # real line parameters
    for j in nps
      S     = hc - randexp()
      L, R  = find_real_int(p, pp, fp, j, S, lhf, w[j], npars)
      p, hc = sample_int(p, pp, fp, j, L, R, S, lhf, npars)
    end

    # hidden factors
    for j in phid
     S     = hc - randexp()
     L, R  = find_nonneg_int(p, pp, fp, fpp, j, S, lhf, w[j], npars)
     p, fp, hc = sample_int(p, pp, fp, fpp, j, L, R, S, lhf, npars)
    end

    #=
    multivariate updates
    =#
    # non hidden factors
    for (i,mvp) in enumerate(mvps)
      S = hc - randexp()
      @views find_rect(p, pp, fp, mvp, nngps[i], S, lhf, w, Lv, Rv, npars)
      p, hc  = sample_rect(p, pp, fp, Lv, Rv, mvp, S, lhf, npars)
    end

    # hidden factors
    for mvp in mvhfs
      S = hc - randexp()
      @views find_rect(p, pp, fp, fpp, mvp, S, lhf, w, Lv, Rv, npars)
      p, fp, hc = sample_rect(p, pp, fp, fpp, Lv, Rv, mvp, S, lhf, npars)
    end

    # log samples
    lthin += 1
    if lthin == nthin
      @inbounds begin
        lit += 1
        setindex!(its,  it, lit)
        setindex!(hlog, hc, lit)
        setindex!(ps,   p,  lit, :)
      end
      lthin = 0
    end

    next!(prog)
  end

  return its, hlog, ps
end






"""
    w_sampler(lhf         ::Function, 
              p           ::Array{Float64,1},
              fp          ::Array{Float64,1},
              nnps        ::Array{Int64,1},
              nps         ::Array{Int64,1},
              phid        ::Array{Int64,1},
              mvps        ::Array{Array{Int64,1},1},
              mvhfs       ::Array{Array{Int64,1},1},
              npars       ::Int64,
              optimal_w   ::Float64,
              screen_print::Int64,
              nburn       ::Int64,
              ntakew      ::Int64,
              winit       ::Float64)

Run slice sampler for burn-in and to estimate appropriate w's.
"""
function w_sampler(lhf         ::Function, 
                   p           ::Array{Float64,1},
                   fp          ::Array{Float64,1},
                   nnps        ::Array{Int64,1},
                   nps         ::Array{Int64,1},
                   phid        ::Array{Int64,1},
                   mvps        ::Array{Array{Int64,1},1},
                   nngps        ::Array{Array{Bool,1},1},
                   mvhfs       ::Array{Array{Int64,1},1},
                   npars       ::Int64,
                   optimal_w   ::Float64,
                   screen_print::Int64,
                   nburn       ::Int64,
                   ntakew      ::Int64,
                   winit       ::Float64)

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

  Lv = Array{Float64,1}(undef, maxmvu)
  Rv = Array{Float64,1}(undef, maxmvu)

  w  = fill(winit, npars)
  ps = Array{Float64,2}(undef, nburn, npars)

  # posterior
  hc = lhf(p, fp)

  # preallocate pp and fpp
  pp  = copy(p)
  fpp = copy(fp)

  prog = Progress(nburn, screen_print, "estimating optimal widths...", 20)

  for it in Base.OneTo(nburn)

    #=
    univariate updates
    =#
    # nonnegative parameters
    for j in nnps
     S     = hc - randexp()
     L, R  = find_nonneg_int(p, pp, fp, j, S, lhf, w[j], npars)
     p, hc = sample_int(p, pp, fp, j, L, R, S, lhf, npars)
    end

    # real line parameters
    for j in nps
      S     = hc - randexp()
      L, R  = find_real_int(p, pp, fp, j, S, lhf, w[j], npars)
      p, hc = sample_int(p, pp, fp, j, L, R, S, lhf, npars)
    end

    # hidden factors
    for j in phid
     S     = hc - randexp()
     L, R  = find_nonneg_int(p, pp, fp, fpp, j, S, lhf, w[j], npars)
     p, fp, hc = sample_int(p, pp, fp, fpp, j, L, R, S, lhf, npars)
    end

    #=
    multivariate updates
    =#
    # non hidden factors
    for (i, mvp) in enumerate(mvps)
      S = hc - randexp()
      @views find_rect(p, pp, fp, mvp, nngps[i], S, lhf, w, Lv, Rv, npars)
      p, hc  = sample_rect(p, pp, fp, Lv, Rv, mvp, S, lhf, npars)
    end

    # hidden factors
    for mvp in mvhfs
      S = hc - randexp()
      @views find_rect(p, pp, fp, fpp, mvp, S, lhf, w, Lv, Rv, npars)
      p, fp, hc = sample_rect(p, pp, fp, fpp, Lv, Rv, mvp, S, lhf, npars)
    end

    @inbounds setindex!(ps, p, it, :)

    next!(prog)
  end

  sps = nburn-ntakew

  ps = ps[(nburn-ntakew+1):nburn,:]

  w = optimal_w .* (reduce(max, ps, dims=1) .- reduce(min, ps, dims=1))
  w = reshape(w, size(w,2))

  return (p, fp, w)::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}}
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
                         postf::Function, 
                         w    ::Float64,
                         npars::Int64)

  unsafe_copyto!(pp, 1, p, 1, npars)

  L::Float64 = pp[j] - w*rand()
  R::Float64 = L + w

  if L <= 0.0
    L = 1e-30
  end

  # left extreme
  pp[j] = L
  while S < postf(pp, fp)
    L -= w
    if L <= 0.0
      L = 1e-30
      break
    end
    pp[j] = L
  end

  # right extreme
  pp[j] = R
  while S < postf(pp, fp)
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
                         postf::Function, 
                         w    ::Float64,
                         npars::Int64)

  unsafe_copyto!(pp,  1, p,  1, npars)
  unsafe_copyto!(fpp, 1, fp, 1, npars)

  L = fpp[j] - w*rand()
  R = L + w

  if L <= 0.0
    L = 1e-30
  end

  # left extreme
  fpp[j] = L
  while S < postf(pp, fpp)
    L -= w
    if L <= 0.0
      L = 1e-30
      break
    end
    fpp[j] = L
  end

  # right extreme
  fpp[j] = R
  while S < postf(pp, fpp)
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
                       postf::Function, 
                       wj   ::Float64,
                       npars::Int64)

  unsafe_copyto!(pp, 1, p, 1, npars)

  L = pp[j] - wj*rand()
  R = L + wj

  # left extreme
  pp[j] = L::Float64
  while S < postf(pp, fp)
    L    -= wj::Float64
    pp[j] = L::Float64
  end

  # right extreme
  pp[j] = R::Float64
  while S < postf(pp, fp)
    R    += wj::Float64
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
                    postf::Function, 
                    npars::Int64)

  @inbounds begin
    unsafe_copyto!(pp, 1, p, 1, npars)

    while true
      pp[j] = (L + rand()*(R-L))::Float64

      hc = postf(pp, fp)::Float64
      if S < hc
        unsafe_copyto!(p, 1, pp, 1, npars)
        return (p, hc)::Tuple{Array{Float64,1}, Float64}
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
                    postf::Function, 
                    npars::Int64)

  @inbounds begin
    unsafe_copyto!(pp,  1, p,  1, npars)
    unsafe_copyto!(fpp, 1, fp, 1, npars)

    while true
      fpp[j] = (L + rand()*(R-L))::Float64

      hc = postf(pp, fpp)::Float64
      if S < hc
        unsafe_copyto!(p,  1, pp,  1, npars)
        unsafe_copyto!(fp, 1, fpp, 1, npars)
        return (p, fp, hc)::Tuple{Array{Float64,1}, Array{Float64,1}, Float64}
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
                   postf::Function, 
                   w    ::Array{Float64,1},
                   Lv   ::Array{Float64,1},
                   Rv   ::Array{Float64,1},
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
      while S < postf(pp, fp)
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
      while S < postf(pp, fp)
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
              mvp  ::Array{Int64,1}, 
              nngp  ::Array{Bool,1},
              S    ::Float64, 
              postf::Function, 
              w    ::Array{Float64,1},
              Lv   ::Array{Float64,1},
              Rv   ::Array{Float64,1},
              npars::Int64)

Returns left and right vertices of an hyper-rectangle for hidden factors.
"""
function find_rect(p    ::Array{Float64,1}, 
                   pp   ::Array{Float64,1}, 
                   fp   ::Array{Float64,1}, 
                   fpp  ::Array{Float64,1},
                   mvp  ::Array{Int64,1}, 
                   S    ::Float64, 
                   postf::Function, 
                   w    ::Array{Float64,1},
                   Lv   ::Array{Float64,1},
                   Rv   ::Array{Float64,1},
                   npars::Int64)

  @inbounds begin

    unsafe_copyto!(pp,  1, p, 1, npars)
    unsafe_copyto!(fpp, 1, fp, 1, npars)

    # randomly start rectangle
    for (i,j) in enumerate(mvp)
      Lv[i] = fpp[j] - w[j]*rand()
      if Lv[i] <= 0.0
        Lv[i] = 1e-30
      end
      Rv[i] = Lv[i] + w[j]
    end

    # left extremes
    for (i,j) in enumerate(mvp)
      fpp[j] = Lv[i]::Float64
      while S < postf(pp, fpp)
        Lv[i] -= w[j]::Float64
        if Lv[i] <= 0.0
          Lv[i] = 1e-30
          break
        end
        fpp[j] = Lv[i]::Float64
      end
    end

    # right extremes
    for (i,j) in enumerate(mvp)
      fpp[j] = Rv[i]::Float64
      while S < postf(pp, fpp)
        Rv[i] += w[j]::Float64
        fpp[j]  = Rv[i]::Float64
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
                     postf::Function, 
                     npars::Int64)

  @inbounds begin

    unsafe_copyto!(pp, 1, p, 1, npars)

    while true

      # multivariate uniform sample
      for (i,j) in enumerate(mvp)
        pp[j] = Lv[i] + rand()*(Rv[i] - Lv[i])
      end

      # posterior
      hc = postf(pp, fp)

      if S < hc
        unsafe_copyto!(p, 1, pp, 1, npars)
        return(p, hc)::Tuple{Array{Float64,1}, Float64}
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
                Lv   ::Array{Float64,1}, 
                Rv   ::Array{Float64,1}, 
                mvp  ::Array{Float64,1},
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
                     S    ::Float64, 
                     postf::Function, 
                     npars::Int64)

  @inbounds begin

    unsafe_copyto!(pp,  1, p,  1, npars)
    unsafe_copyto!(fpp, 1, fp, 1, npars)

    while true

      # multivariate uniform sample
      for (i,j) in enumerate(mvp)
        fpp[j] = Lv[i] + rand()*(Rv[i] - Lv[i])
      end

      # posterior
      hc = postf(pp, fpp)

      if S < hc
        unsafe_copyto!(p, 1, pp, 1, npars)
        unsafe_copyto!(fp, 1, fpp, 1, npars)
        return(p, fp, hc)::Tuple{Array{Float64,1}, Array{Float64,1}, Float64}
      end

      for (i,j) in enumerate(mvp)
        if fpp[j] < fp[j]
            Lv[i] = fpp[j]::Float64
        else
            Rv[i] = fpp[j]::Float64
        end
      end

    end
  end

  return nothing
end

