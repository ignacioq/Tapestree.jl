#=

MCMC Metropolis-Hastings sampling

Ignacio Quintero Mächler

t(-_-t)

29 04 2021

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
    mcmcmh(lhf         ::Function, 
           p           ::DArray{Array{Float64,1},1,Array{Array{Float64,1},1}},
           fp          ::DArray{Array{Float64,1},1,Array{Array{Float64,1},1}},
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

Run slice-sampling Markov Chain given posterior function.
"""
function mcmcmh(lhf         ::Function, 
                p           ::DArray{Array{Float64,1},1,Array{Array{Float64,1},1}},
                fp          ::DArray{Array{Float64,1},1,Array{Array{Float64,1},1}},
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

  # temperature
  t, o = make_temperature(0.0, ncch)

  # make SharedArray for temperature order `o`
  o = SharedArray(o)

  # run burnin phase
  p, fp, tn, o, t = 
    mcmcmh_burn(lhf, p, fp, nnps, nps, phid, nburn, ncch, nswap, 
      o, t, tni, tune_int, npars, screen_print, obj_ar)

  temperature!(t, o, dt)

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
  I = Float64.(il[:,1])

  d_closeall()

  # save samples
  R = hcat(I, H, P)

  return R
end





