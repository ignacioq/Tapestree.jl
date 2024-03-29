#=

MCMC Metropolis-Hastings utilities

Ignacio Quintero Mächler

t(-_-t)

23 04 2021

=#




"""
    mcmcmh_burn(lhf         ::Function, 
                p           ::Array{Array{Float64,1},1},
                fp          ::Array{Array{Float64,1},1},
                nnps        ::Array{Int64,1},
                nps         ::Array{Int64,1},
                phid        ::Array{Int64,1},
                nburn       ::Int64,
                ncch        ::Int64,
                nswap       ::Int64,
                dt          ::Float64
                tni         ::Float64,
                tune_int    ::Int64,
                npars       ::Int64,
                screen_print::Int64)

Run Metropolis-Hastings burn-in phase.
"""
function mcmcmh_burn(lhf         ::Function, 
                     p           ::Array{Array{Float64,1},1},
                     fp          ::Array{Array{Float64,1},1},
                     nnps        ::Array{Int64,1},
                     nps         ::Array{Int64,1},
                     phid        ::Array{Int64,1},
                     nburn       ::Int64,
                     ncch        ::Int64,
                     nswap       ::Int64,
                     dt          ::Float64,
                     tni         ::Float64,
                     tune_int    ::Int64,
                     npars       ::Int64,
                     screen_print::Int64,
                     obj_ar      ::Float64)

  # temperature
  t, o = make_temperature(dt, ncch)

  # preallocate pp and fpp and tuning
  pp  = copy(p[1])
  fpp = copy(fp[1])
  lfp = lastindex(fpp)
  tn  = fill(tni, npars)

  # make Array for current posteriors
  lhc  = [lhf(p[c], fp[c], t[o[c]]) for c in Base.OneTo(ncch)]

  # make array of arrays for tuning and acceptance rates
  ltn = 0
  lup = 0.0
  lac = [zeros(Float64, npars) for i in Base.OneTo(ncch)]

  scalef = makescalef(obj_ar)

  lswap = 0

  prog = Progress(nburn, screen_print, "mcmc-mh burn-in...", 20)

  # start burn-in iterations
  for it in Base.OneTo(nburn) 

    # updates
    for c in Base.OneTo(ncch)
      lhc[c] = par_cycle(lhf, lhc[c], p[c], fp[c], pp, fpp, lac[c], 
                 nnps, nps, phid, tn, t[o[c]], npars, lfp)
    end

    # log tuning parameters only for T = 1
    lup += 1.0
    ltn += 1
    if ltn == tune_int
      lac1 = lac[findfirst(isone, o)]
      for j in nnps
        tn[j] = scalef(tn[j], lac1[j]/lup)
      end
      for j in nps
        tn[j] = scalef(tn[j], lac1[j]/lup)
      end
      for j in phid
        tn[j] = scalef(tn[j], lac1[j]/lup)
      end
      ltn = 0
    end

    # swap chains
    if ncch > 1
      lswap += 1
      if lswap == nswap
        swap_chains!(o, t, lhc)
        lac = lac[o]
        lswap = 0
      end
    end

    next!(prog)
  end

  return p, fp, tn, o, t
end






"""
    mcmcmh_burn(lhf         ::Function, 
                p           ::DArray{Array{Float64,1},1,Array{Array{Float64,1},1}},
                fp          ::DArray{Array{Float64,1},1,Array{Array{Float64,1},1}},
                nnps        ::Array{Int64,1},
                nps         ::Array{Int64,1},
                phid        ::Array{Int64,1},
                nburn       ::Int64,
                ncch        ::Int64,
                nswap       ::Int64,
                dt          ::Float64,
                tni         ::Float64,
                tune_int    ::Int64,
                npars       ::Int64,
                screen_print::Int64,
                obj_ar      ::Float64)

Run parallel Metropolis-Hastings burn-in phase.
"""
function mcmcmh_burn(lhf         ::Function, 
                     p           ::DArray{Array{Float64,1},1,Array{Array{Float64,1},1}},
                     fp          ::DArray{Array{Float64,1},1,Array{Array{Float64,1},1}},
                     nnps        ::Array{Int64,1},
                     nps         ::Array{Int64,1},
                     phid        ::Array{Int64,1},
                     nburn       ::Int64,
                     ncch        ::Int64,
                     nswap       ::Int64,
                     o           ::SharedArray{Int64,1}, 
                     t           ::Array{Float64,1},
                     tni         ::Float64,
                     tune_int    ::Int64,
                     npars       ::Int64,
                     screen_print::Int64,
                     obj_ar      ::Float64)

  # Make distributed array for iterations per chain `ipc`
  ipc = SharedArray{Int64,1}(zeros(Int64,ncch))

  # make Array for current posteriors
  lhc  = [lhf(p[c], fp[c], t[o[c]]) for c in Base.OneTo(ncch)]
  lhc = SharedArray(lhc)

  # pre-estimate chain swaps
  tns   = fld(nburn, nswap)
  allsw = [wchains(o) for i in Base.OneTo(tns)]
  swr   = SharedArray{Bool,1}(fill(false, tns))
  swl   = SharedArray{Bool,1}(fill(false, tns))

  # tuning array
  tn = [fill(tni, npars) for i in Base.OneTo(ncch)]
  tn = distribute(tn)

  # make array of arrays for tuning and acceptance rates
  lup = zeros(Float64, ncch)
  lup = SharedArray(lup)
  lac = [zeros(Float64, npars) for i in Base.OneTo(ncch)]
  lac = distribute(lac)

  scalef = makescalef(obj_ar)

  lfp = lastindex(fp[1])

  # start parallel mc3
  @sync @distributed for c in Base.OneTo(ncch)

    # swap variables
    lswap, ws = 0, 0
    prog = Progress(nburn, screen_print, "mcmc-mh burn-in...", 20)

    # preallocate pp and fpp and tuning
    pp  = copy(p[c])
    fpp = copy(fp[c])

    ltn = 0

    # start burn-in iterations
    for it in Base.OneTo(nburn) 

      lhc[c] = par_cycle(lhf, lhc[c], p[c], fp[c], pp, fpp, lac[c], 
                 nnps, nps, phid, tn[c], t[o[c]], npars, lfp)

      # log tuning parameters only for T = 1
      lup[c] += 1.0
      ltn += 1
      if ltn == tune_int
        tnc  = tn[c]
        lacc = lac[c]
        for j in nnps
          tnc[j] = scalef(tnc[j], lacc[j]/lup[c])
        end
        for j in nps
          tnc[j] = scalef(tnc[j], lacc[j]/lup[c])
        end
        for j in phid
          tnc[j] = scalef(tnc[j], lacc[j]/lup[c])
        end
        ltn = 0
      end

      # log chain current iter
      ipc[c] = it

      ## swap chains
      lswap += 1
      if lswap == nswap
        ws += 1
        j, k = allsw[ws]
        cij = c == j   # c is j
        cik = c == k   # c is k
        # swapper
        if cij
          while true
            ipc[c] == ipc[k] && swl[ws] && break
          end
          swr[ws] = swap_chains!(j, k, o, lhc)
        end
        # swappee
        if cik
          swl[ws] = true
          while true
            ipc[c] == ipc[j] && swr[ws] && break
          end
        end
        lswap = 0
      end
      next!(prog)
    end
  end

  return p, fp, tn[1], o, t
end




"""
    mcmcmh_mcmc(lhf         ::Function, 
                p           ::Array{Array{Float64,1},1},
                fp          ::Array{Array{Float64,1},1},
                nnps        ::Array{Int64,1},
                nps         ::Array{Int64,1},
                phid        ::Array{Int64,1},
                niter       ::Int64,
                nthin       ::Int64,
                ncch        ::Int64,
                nswap       ::Int64,
                tn          ::Array{Float64,1},
                o           ::Array{Int64,1}, 
                t           ::Array{Float64,1},
                npars       ::Int64,
                out_file    ::String,
                pardic      ::Dict{String, Int64},
                screen_print::Int64)

Run Metropolis-Hastings burn-in phase.
"""
function mcmcmh_mcmc(lhf         ::Function, 
                     p           ::Array{Array{Float64,1},1},
                     fp          ::Array{Array{Float64,1},1},
                     nnps        ::Array{Int64,1},
                     nps         ::Array{Int64,1},
                     phid        ::Array{Int64,1},
                     niter       ::Int64,
                     nthin       ::Int64,
                     ncch        ::Int64,
                     nswap       ::Int64,
                     tn          ::Array{Float64,1},
                     o           ::Array{Int64,1}, 
                     t           ::Array{Float64,1},
                     npars       ::Int64,
                     out_file    ::String,
                     pardic      ::Dict{String, Int64},
                     screen_print::Int64)

  nlogs = fld(niter, nthin)

  #preallocate logging arrays
  il =  Array{Float64,1}(undef, nlogs)
  hl =  Array{Float64,2}(undef, nlogs, ncch)
  pl = [Array{Float64,2}(undef, nlogs, npars) for i in Base.OneTo(ncch)]

  # preallocate chains order
  ol   = Array{Int64,2}(undef, nlogs, ncch)

  lthin, lit, lswap = 0, 0, 0

  # preallocate pp and fpp and tuning
  pp  = copy(p[1])
  fpp = copy(fp[1])
  lfp = lastindex(fpp)

  # starting posteriors
  lhc  = [lhf(p[c], fp[c], t[o[c]]) for c in Base.OneTo(ncch)]

  prog = Progress(niter, screen_print, "running mcmc-mh...", 20)

  # make header of file
  pdc = sort!(collect(pardic), by = x -> x[2])
  hh  = "Iteration\tPosterior"
  for (k,v) in pdc
    hh *= "\t"*k
  end
  hh *= "\n"

  open(out_file*"_flush.log", "w") do of

    write(of, hh)
    flush(of)

    # start burn-in iterations
    for it in Base.OneTo(niter) 

      # updates
      for c in Base.OneTo(ncch)
        lhc[c] = par_cycle(lhf, lhc[c], p[c], fp[c], pp, fpp, 
                   nnps, nps, phid, tn, t[o[c]], npars, lfp)
      end

      # log samples
      lthin += 1
      if lthin === nthin
        # write to table
        @inbounds begin
          lit += 1
          setindex!(il,  it,  lit)
          setindex!(hl, lhc, lit, :)
          for c in Base.OneTo(ncch)
            setindex!(pl[c], p[c], lit, :)
          end
          setindex!(ol, o, lit, :)
        end

        # write to file (only chain 1)
        o1 = findfirst(isone, o)
        p1 = p[o1]
        sts = string(it,"\t",lhc[o1],"\t")
        for i in Base.OneTo(npars)
          sts *= string(p1[i], "\t")
        end
        sts *= "\n"
        write(of, sts)
        flush(of)

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
  end

  return il, hl, pl, ol
end





"""
    mcmcmh_mcmc(lhf         ::Function, 
                p           ::DArray{Array{Float64,1},1,Array{Array{Float64,1},1}},
                fp          ::DArray{Array{Float64,1},1,Array{Array{Float64,1},1}},
                nnps        ::Array{Int64,1},
                nps         ::Array{Int64,1},
                phid        ::Array{Int64,1},
                niter       ::Int64,
                nthin       ::Int64,
                ncch        ::Int64,
                nswap       ::Int64,
                tn          ::Array{Float64,1},
                o           ::SharedArray{Int64,1}, 
                t           ::Array{Float64,1},
                npars       ::Int64,
                out_file    ::String,
                pardic      ::Dict{String, Int64},
                screen_print::Int64)

Run Metropolis-Hastings burn-in phase.
"""
function mcmcmh_mcmc(lhf         ::Function, 
                     p           ::DArray{Array{Float64,1},1,Array{Array{Float64,1},1}},
                     fp          ::DArray{Array{Float64,1},1,Array{Array{Float64,1},1}},
                     nnps        ::Array{Int64,1},
                     nps         ::Array{Int64,1},
                     phid        ::Array{Int64,1},
                     niter       ::Int64,
                     nthin       ::Int64,
                     ncch        ::Int64,
                     nswap       ::Int64,
                     tn          ::Array{Float64,1},
                     o           ::SharedArray{Int64,1}, 
                     t           ::Array{Float64,1},
                     npars       ::Int64,
                     out_file    ::String,
                     pardic      ::Dict{String, Int64},
                     screen_print::Int64)

  nlogs = fld(niter, nthin)

  #preallocate logging arrays
  il =  SharedArray{Int64,2}(zeros(Int64,nlogs, ncch))
  hl =  SharedArray{Float64,2}(zeros(nlogs, ncch))
  pl = [Array{Float64,2}(undef, nlogs, npars) for i in Base.OneTo(ncch)]
  pl = distribute(pl)

  # Make distributed array for iterations per chain `ipc`
  ipc = SharedArray{Int64,1}(zeros(Int64,ncch))

  # preallocate chains order
  ol   = SharedArray{Int64,2}(zeros(Int64, nlogs, ncch))

  # length fp
  lfp = length(fp[1])

  # starting posteriors
  lhc = SharedArray([lhf(p[c], fp[c], t[o[c]]) for c in Base.OneTo(ncch)])

  # pre-estimate chain swaps 
  tns   = fld(niter, nswap)
  allsw = [wchains(o) for i in Base.OneTo(tns)]
  swr   = SharedArray{Bool,1}(fill(false, tns))
  swl   = SharedArray{Bool,1}(fill(false, tns))

  # make file header
  pdc = sort!(collect(pardic), by = x -> x[2])
  hh  = "Iteration\tPosterior"
  for (k,v) in pdc
    hh *= "\t"*k
  end
  hh *= "\tchain_t\n"

  # start parallel slice-sampling
  @sync @distributed for c in Base.OneTo(ncch)

    open(string(out_file,"_",c,"_flush.log"), "w") do of

      write(of, hh)
      flush(of)

      # log variables
      lthin, lit, lswap, ws = 0, 0, 0, 0
      prog = Progress(niter, screen_print, "running mcmc-mh...", 20)

      # preallocate pp and fpp
      pp  = copy(p[c])
      fpp = copy(fp[c])

     # start burn-in iterations
      for it in Base.OneTo(niter) 

        # updates
        lhc[c] = par_cycle(lhf, lhc[c], p[c], fp[c], pp, fpp, 
                   nnps, nps, phid, tn, t[o[c]], npars, lfp)

        # log samples
        lthin += 1
        if lthin === nthin
          # write to table
          @inbounds begin
            lit += 1
            il[lit, c]    = it
            hl[lit, c]    = lhc[c]
            ol[lit, c]    = o[c]
            pl[c][lit, :] = p[c]
          end

          # write to file
          sts = string(it,"\t",lhc[c],"\t")
          for i in Base.OneTo(npars)
            sts *= string(p[c][i], "\t")
          end
          sts *= string(o[c], "\n")
          write(of, sts)
          flush(of)

          lthin = 0
        end

        # log chain current iter
        ipc[c] = it

        ## swap chains
        lswap += 1
        if lswap == nswap
          ws += 1
          j, k = allsw[ws]
          cij = c == j   # c is j
          cik = c == k   # c is k
          # swapper
          if cij
            while true
              ipc[c] == ipc[k] && swl[ws] && break
            end
            swr[ws] = swap_chains!(j, k, o, t, lhc)
          end
          # swappee
          if cik
            swl[ws] = true
            while true
              ipc[c] == ipc[j] && swr[ws] && break
            end
          end
          lswap = 0
        end
        next!(prog)
      end
    end
  end

  return il, hl, pl, ol
end





"""
    mulupt(p::Float64, tn::Float64)

Multiplicative parameter window move.
"""
mulupt(p::Float64, tn::Float64) = p * exp((rand() - 0.5) * tn)




"""
    addupt(p::Float64, tn::Float64)

Gaussian parameter window move.
"""
addupt(p::Float64, tn::Float64) = p + randn() * tn




"""
    par_cycle(lhf  ::Function, 
              lhc  ::Float64,
              p    ::Array{Float64,1},
              fp   ::Array{Float64,1},
              pp   ::Array{Float64,1}, 
              fpp  ::Array{Float64,1}, 
              lac  ::Array{Float64,1},
              nnps ::Array{Int64,1},
              nps  ::Array{Int64,1},
              phid ::Array{Int64,1},
              tn   ::Array{Float64,1},
              ti   ::Float64)

A full parameter update MH cycle.
"""
function par_cycle(lhf  ::Function, 
                   lhc  ::Float64,
                   p    ::Array{Float64,1},
                   fp   ::Array{Float64,1},
                   pp   ::Array{Float64,1}, 
                   fpp  ::Array{Float64,1}, 
                   lac  ::Array{Float64,1},
                   nnps ::Array{Int64,1},
                   nps  ::Array{Int64,1},
                   phid ::Array{Int64,1},
                   tn   ::Array{Float64,1},
                   ti   ::Float64,
                   npars::Int64,
                   lfp  ::Int64)

  # nonnegative parameters
  for j in shuffle!(nnps)
    unsafe_copyto!(pp, 1, p, 1, npars)

    pp[j] = mulupt(p[j], tn[j])
    lhp   = lhf(pp, fp, ti)

    if -randexp() < (lhp - lhc[1] + log(pp[j]/p[j]))
      lhc = lhp
      unsafe_copyto!(p, 1, pp, 1, npars)

      lac[j] += 1.0
    end
  end

  # real line parameters
  for j in shuffle!(nps)
    unsafe_copyto!(pp, 1, p, 1, npars)

    pp[j] = addupt(p[j], tn[j])
    lhp   = lhf(pp, fp, ti)

    if -randexp() < (lhp - lhc)
      lhc = lhp
      unsafe_copyto!(p, 1, pp, 1, npars)

      lac[j] += 1.0
    end
  end

  # hidden factors
  for j in shuffle!(phid)

    unsafe_copyto!(fpp, 1, fp, 1, lfp)
    unsafe_copyto!(pp, 1, p, 1, npars)

    fpp[j] = mulupt(fp[j], tn[j])
    lhp    = lhf(p, fpp, ti)

    if -randexp() < (lhp - lhc + log(fpp[j]/fp[j]))
      lhc = lhp
      unsafe_copyto!(fp, 1, fpp, 1, lfp)
      unsafe_copyto!(p, 1, pp, 1, npars)
      lac[j] += 1.0
    end
  end

  return lhc
end





"""
    par_cycle(lhf  ::Function, 
              lhc  ::Float64,
              p    ::Array{Float64,1},
              fp   ::Array{Float64,1},
              pp   ::Array{Float64,1}, 
              fpp  ::Array{Float64,1}, 
              nnps ::Array{Int64,1},
              nps  ::Array{Int64,1},
              phid ::Array{Int64,1},
              tn   ::Array{Float64,1},
              ti   ::Float64)

A full parameter update MH cycle.
"""
function par_cycle(lhf  ::Function, 
                   lhc  ::Float64,
                   p    ::Array{Float64,1},
                   fp   ::Array{Float64,1},
                   pp   ::Array{Float64,1}, 
                   fpp  ::Array{Float64,1}, 
                   nnps ::Array{Int64,1},
                   nps  ::Array{Int64,1},
                   phid ::Array{Int64,1},
                   tn   ::Array{Float64,1},
                   ti   ::Float64,
                   npars::Int64,
                   lfp  ::Int64)


  # nonnegative parameters
  for j in shuffle!(nnps)
    unsafe_copyto!(pp, 1, p, 1, npars)

    pp[j] = mulupt(p[j], tn[j])
    lhp   = lhf(pp, fp, ti)

    if -randexp() < (lhp - lhc[1] + log(pp[j]/p[j]))
      lhc     = lhp
      unsafe_copyto!(p, 1, pp, 1, npars)
    end
  end

  # real line parameters
  for j in shuffle!(nps)
    unsafe_copyto!(pp, 1, p, 1, npars)

    pp[j] = addupt(p[j], tn[j])
    lhp   = lhf(pp, fp, ti)

    if -randexp() < (lhp - lhc)
      lhc    = lhp
      unsafe_copyto!(p, 1, pp, 1, npars)
    end
  end

  # hidden factors
  for j in shuffle!(phid)
    unsafe_copyto!(fpp, 1, fp, 1, lfp)
    unsafe_copyto!(pp, 1, p, 1, npars)

    fpp[j] = mulupt(fp[j], tn[j])
    lhp    = lhf(p, fpp, ti)

    if -randexp() < (lhp - lhc + log(fpp[j]/fp[j]))
      lhc    = lhp
      unsafe_copyto!(fp, 1, fpp, 1, lfp)
      unsafe_copyto!(p, 1, pp, 1, npars)
    end
  end

  return lhc
end




"""
    makescalef(obj_ar::Float64)

Make scaling function given the objective acceptance rates.
"""
function makescalef(obj_ar::Float64)
  inar::Float64 = 1.0/(1.0 - obj_ar)
  iobj_ar::Float64 = 1.0/obj_ar

  f(window::Float64, rate::Float64) = let obj_ar = obj_ar, inar = inar, iobj_ar = iobj_ar
    if rate > obj_ar
      window *= (1.0 + (rate - obj_ar)*inar)::Float64
    else
      window /= (2.0 - rate*iobj_ar)::Float64
    end
    
    return window::Float64
  end

  return f
end



