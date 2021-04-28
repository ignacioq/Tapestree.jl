#=

Slice sampling for EGeoHiSSE

Ignacio Quintero Mächler

t(-_-t)

November 20 2017

=#





"""
    mcmcmh(lhf         ::Function, 
           p           ::Array{Array{Float64,1},1},
           fp          ::Array{Array{Float64,1},1},
           nnps        ::Array{Int64,1},
           nps         ::Array{Int64,1},
           phid        ::Array{Int64,1},
           npars       ::Int64,
           niter       ::Int64,
           nthin       ::Int64,
           nburn       ::Int64,
           nswap       ::Int64,
           ncch        ::Int64,
           tni         ::Float64,
           tune_int    ::Int64,
           dt          ::Float64,
           screen_print::Int64)

Run slice-sampling Markov Chain given posterior function.
"""
function mcmcmh(lhf         ::Function, 
                p           ::Array{Array{Float64,1},1},
                fp          ::Array{Array{Float64,1},1},
                nnps        ::Array{Int64,1},
                nps         ::Array{Int64,1},
                phid        ::Array{Int64,1},
                npars       ::Int64,
                niter       ::Int64,
                nthin       ::Int64,
                nburn       ::Int64,
                nswap       ::Int64,
                ncch        ::Int64,
                tni         ::Float64,
                tune_int    ::Int64,
                dt          ::Float64,
                screen_print::Int64,
                obj_ar      ::Float64)

  # run burnin phase
  p, fp, tn, o, t = 
    mcmcmh_burn(lhf, p, fp, nnps, nps, phid, nburn, ncch, nswap, 
      dt, tni, tune_int, npars, screen_print, obj_ar)

  # run mcmc
  il, hl, pl, ol = 
    mcmcmh_mcmc(lhf, p, fp, nnps, nps, phid, niter, nthin, ncch, nswap, 
      tn, o, t, npars, screen_print)

  # choose cold chain (is equal to 1)
  P = Array{Float64,2}(undef, size(hl,1), npars)
  H = Array{Float64,1}(undef, size(hl,1))
  @inbounds begin
    for i in Base.OneTo(size(hl,1))
      @views ii = findfirst(x -> isone(x), ol[i,:])
      P[i,:] = pl[ii][i,:]
      H[i]   = hl[i,ii]
    end
  end

  # save samples
  R = hcat(il, H, P)

  return R
end






"""
    slice_sampler(tip_val     ::Dict{Int64,Array{Float64,1}},
                  ed          ::Array{Int64,2},
                  el          ::Array{Float64,1},
                  x           ::Array{Float64,1},
                  y           ::Array{Float64},
                  cov_mod     ::String,
                  out_file    ::String,
                  h           ::Int64;
                  constraints ::NTuple{N,String}  = (" ",),
                  niter       ::Int64             = 10_000,
                  nthin       ::Int64             = 10,
                  λpriors     ::Float64           = .1,
                  μpriors     ::Float64           = .1,
                  gpriors     ::Float64           = .1,
                  lpriors     ::Float64           = .1,
                  qpriors     ::Float64           = .1,
                  βpriors     ::NTuple{2,Float64} = (0.0, 5.0),
                  optimal_w   ::Float64           = 0.8,
                  screen_print::Int64             = 5) where {N}

Parallel run slice-sampling Markov Chain given posterior function.
"""
function slice_sampler(lhf         ::Function, 
                       p           ::DArray{Array{Float64,1},1,Array{Array{Float64,1},1}},
                       fp          ::DArray{Array{Float64,1},1,Array{Array{Float64,1},1}},
                       nnps        ::Array{Int64,1},
                       nps         ::Array{Int64,1},
                       phid        ::Array{Int64,1},
                       mvps        ::Array{Array{Int64,1},1},
                       nngps       ::Array{Array{Bool,1},1},
                       mvhfs       ::Array{Array{Int64,1},1},
                       hfgps       ::Array{Array{Bool,1},1},
                       npars       ::Int64,
                       niter       ::Int64,
                       nthin       ::Int64,
                       nburn       ::Int64,
                       ntakew      ::Int64,
                       nswap       ::Int64,
                       ncch        ::Int64,
                       winit       ::Float64,
                       optimal_w   ::Float64,
                       dt          ::Float64,
                       screen_print::Int64)

  # estimate optimal w
  p, fp, w, o, t = 
    w_sampler(lhf, p, fp, nnps, nps, phid, mvps, nngps, mvhfs, hfgps,
      npars, optimal_w, screen_print, nburn, ntakew, nswap, ncch, winit, dt)

  # slice-sampler
  il, hl, pl, ol = 
    loop_slice_sampler(lhf, p, fp, nnps, nps, phid, mvps, nngps, mvhfs, hfgps, 
      w, npars, niter, nthin, nswap, ncch, o, t, screen_print)

  # choose cold chain (is equal to 1)
  P = Array{Float64,2}(undef, size(hl,1), npars)
  H = Array{Float64,1}(undef, size(hl,1))
  @inbounds begin
    for i in Base.OneTo(size(hl,1))
      @views ii = findfirst(x -> isone(x), ol[i,:])
      P[i,:] = pl[ii][i,:]
      H[i]   = hl[i,ii]
    end
  end
  I = Float64.(il[:,1])

  d_closeall()

  # save samples
  R = hcat(I, H, P)

  return R
end


