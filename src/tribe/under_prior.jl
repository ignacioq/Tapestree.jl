#=

MCMC for biogeographic competition model

Ignacio Quintero Mächler

t(-_-t)

April 27 2017

=#





"""
    tribe_mcmc(...)

Run MCMC for trait and range interspecitific biogeographic evolution (tribe)
**under the prior**.
"""
function tribe_mcmc(out_file::String;
                    niter   ::Int64             = 50_000,
                    nthin   ::Int64             = 100,
                    nburn   ::Int64             = 50_000,
                    ωxprior ::NTuple{2,Float64} = (0.,10.),
                    ω1prior ::NTuple{2,Float64} = (1.,1.),
                    ω0prior ::Float64           = 1e-1,
                    σ²prior ::Float64           = 1e-1,
                    λprior  ::Float64           = 1e-1,
                    weight  ::NTuple{4,Float64} = (0.15,0.05,0.02,0.02),
                    σ²i     ::Float64           = 1.,
                    ωxi     ::Float64           = 0.,
                    ω1i     ::Float64           = 0.01,
                    ω0i     ::Float64           = 0.01,
                    λ1i     ::Float64           = 1.0,
                    λ0i     ::Float64           = 0.2,
                    fix_ωx  ::Bool              = false,
                    fix_ω1  ::Bool              = false,
                    fix_ω0  ::Bool              = false)

  parvec = [1:6...]
  np = length(parvec)

  # add to parameter update vector according to weights
  append!(parvec, fill(1,ceil(Int64,np*weight[1])))
  append!(parvec, fill(2,ceil(Int64,np*weight[2])))
  append!(parvec, repeat(3:4,  inner = ceil(Int64,np*weight[3])))
  append!(parvec, repeat(5:6,  inner = ceil(Int64,np*weight[4])))

  # if fixed ωx, ω1 or ω0, then remove
  fix_ωx && filter!(x -> x ≠ 2, parvec)
  fix_ω1 && filter!(x -> x ≠ 3, parvec)
  fix_ω0 && filter!(x -> x ≠ 4, parvec)

  # burning phase
  llc, prc, σ²c, ωxc, ω1c, ω0c, λ1c, λ0c, ptn = 
    burn_tribe(σ²i, ωxi, ω1i, ω0i, λ1i, λ0i, 
               λprior, ωxprior, ω1prior, ω0prior, σ²prior, 
               np, parvec, nburn)

  # log for nthin
  lthin = 0

  # progress bar
  p = Progress(niter, dt=5, desc="mcmc...", barlen=20, color=:green)

  # write to file
  open(out_file*".log", "w")

  open(out_file*".log", "a") do f

    print(f, "iteration\tlikelihood\tprior\tomega_x\tomega_1\tomega_0\tsigma2\tlambda_1\tlambda_0\tcollision_probability\n")

    #=
    start MCMC
    =#

    for it = Base.OneTo(niter)

      # update vector order
      # Update vector
      shuffle!(parvec)

      for up = parvec

        # update λ1
        if up == 5

          # update λ
          λ1p = mulupt(λ1c, ptn[5])::Float64

          llr = 0.0

          prr = llrdexp_x(λ1p, λ1c, λprior)

          if -Random.randexp()()() < (llr + prr + log(λ1p/λ1c))
            llc    += llr::Float64
            prc    += prr::Float64
            λ1c     = λ1p::Float64
          end

        # update λ0
        elseif up == 6

          # update λ
          λ0p = mulupt(λ0c, ptn[6])::Float64

          llr = 0.0

          prr = llrdexp_x(λ0p, λ0c, λprior)

          if -Random.randexp()()() < (llr + prr + log(λ0p/λ0c))
            llc     += llr::Float64
            prc     += prr::Float64
            λ0c      = λ0p::Float64
          end

        # if σ² is updated
        elseif up == 1

          σ²p = mulupt(σ²c, ptn[1])::Float64

          #likelihood ratio
          if ωxc >= 0.0
            llr = 0.0
          else
            llr = 0.0
          end
          
          # prior ratio
          prr = llrdexp_x(σ²p, σ²c, σ²prior)

          if -Random.randexp()()() < (llr + prr + log(σ²p/σ²c))
            llc += llr::Float64
            prc += prr::Float64
            σ²c  = σ²p::Float64
          end

        #update ωx
        elseif up == 2

          ωxp = addupt(ωxc, ptn[2])::Float64

          #likelihood ratio

          llr = 0.0

          # prior ratio
          prr = llrdnorm_x(ωxp, ωxc, ωxprior[1], ωxprior[2])

          if -Random.randexp()()() < (llr + prr)
            llc += llr::Float64
            prc += prr::Float64
            ωxc  = ωxp::Float64
          end

        #update ω1
        elseif up == 3

          ω1p = addupt(ω1c, ptn[3])::Float64

          # proposal likelihood and prior
          llr = 0.0

          # prior ratio
          prr = llrdbeta_x(ω1p, ω1c, ω1prior[1], ω1prior[2])

          if -Random.randexp()()() < (llr + prr)
            llc += llr::Float64
            prc += prr::Float64
            ω1c  = ω1p::Float64
          end

        # update ω0
        else

          ω0p = mulupt(ω0c, ptn[4])

          llr = 0.0

          # prior ratio
          prr = llrdexp_x(ω0p, ω0c, ω0prior)

          if -Random.randexp()()() < (llr + prr + log(ω0p/ω0c))
            llc += llr
            prc += prr
            ω0c  = ω0p
          end
        end

      end # end parameter loop

      # log parameters
      lthin += 1
      if lthin == nthin
        pci = Pc(f_λ(λ1c,ω1c,1.0), f_λ(λ0c,ω0c,1.0), 0.01)
        print(f, it,"\t", llc, "\t", prc,"\t",ωxc,"\t",ω1c,"\t",ω0c,"\t",
             σ²c,"\t",λ1c,"\t",λ0c,"\t",pci,"\n")
        lthin = 0
      end

      next!(p)
    end # end MCMC
  end # end print loop

  return llc, prc, ωxc, ω1c, ω0c, σ²c, λ1c, λ0c
end





"""
    burn_compete(...)

Burning & adaptive phase for MCMC **under the prior**.
"""
function burn_tribe(σ²c     ::Float64,
                    ωxc     ::Float64,
                    ω1c     ::Float64,
                    ω0c     ::Float64,
                    λ1c     ::Float64,
                    λ0c     ::Float64,
                    λprior  ::Float64,
                    ωxprior ::NTuple{2,Float64},
                    ω1prior ::NTuple{2,Float64},
                    ω0prior ::Float64,
                    σ²prior ::Float64,
                    np      ::Int64,
                    parvec  ::Array{Int64,1},
                    nburn   ::Int64;
                    obj_ar  ::Float64 = 0.234,
                    tune_int::Int64   = 100)

  llc = 0.0

  prc = allλpr(  λ1c, λ0c, λprior)            +
        logdexp( σ²c, σ²prior)                +
        logdnorm(ωxc, ωxprior[1], ωxprior[2]) +
        logdbeta(ω1c, ω1prior[1], ω1prior[2]) +
        logdexp( ω0c, ω0prior)

  # make scaling function
  scalef = makescalef(obj_ar)

  # rest of tuning parameters
  ptn = fill(.1, np) 

  # initialize acceptance log
  ltn = zeros(Int64, np)
  lup = zeros(Int64, np)
  lac = zeros(Int64, np)

  # progress bar
  p = Progress(nburn, dt=5, desc="burn...", barlen=20, color=:green)

  # print number of parameters
  print_with_color(:green,
    "\nωx updates per iter = ", lastindex(filter(x -> x == 2,parvec)),
    "\nω1 updates per iter = ", lastindex(filter(x -> x == 3,parvec)),
    "\nω0 updates per iter = ", lastindex(filter(x -> x == 4,parvec)),
    "\nσ² updates per iter = ", lastindex(filter(x -> x == 1,parvec)),
    "\nλ1 updates per iter = ", lastindex(filter(x -> x == 5,parvec)),
    "\nλ0 updates per iter = ", lastindex(filter(x -> x == 6,parvec)),
    "\nθ  updates per iter = ", length(parvec), "\n")

  #start burnin
  for it = Base.OneTo(nburn)

    # Update vector
    shuffle!(parvec)

    for up = parvec

      # update λ1
      if up == 5

        # update λ
        λ1p = mulupt(λ1c, ptn[5])::Float64

        llr = 0.0

        prr = llrdexp_x(λ1p, λ1c, λprior)

        if -Random.randexp()()() < (llr + prr + log(λ1p/λ1c))
          llc    += llr::Float64
          prc    += prr::Float64
          λ1c     = λ1p::Float64
          lac[5] += 1
        end

      # update λ0
      elseif up == 6

        # update λ
        λ0p = mulupt(λ0c, ptn[6])::Float64

        llr = 0.0

        prr = llrdexp_x(λ0p, λ0c, λprior)

        if -Random.randexp()()() < (llr + prr + log(λ0p/λ0c))
          llc     += llr::Float64
          prc     += prr::Float64
          λ0c      = λ0p::Float64
          lac[6] += 1
        end

      # if σ² is updated
      elseif up == 1

        σ²p = mulupt(σ²c, ptn[1])::Float64

        #likelihood ratio
        if ωxc >= 0.0
          llr = 0.0
        else
          llr = 0.0
        end
        
        # prior ratio
        prr = llrdexp_x(σ²p, σ²c, σ²prior)

        if -Random.randexp()()() < (llr + prr + log(σ²p/σ²c))
          llc += llr::Float64
          prc += prr::Float64
          σ²c  = σ²p::Float64
          lac[1] += 1
        end

      #update ωx
      elseif up == 2

        ωxp = addupt(ωxc, ptn[2])::Float64

        #likelihood ratio

        llr = 0.0

        # prior ratio
        prr = llrdnorm_x(ωxp, ωxc, ωxprior[1], ωxprior[2])

        if -Random.randexp()()() < (llr + prr)
          llc += llr::Float64
          prc += prr::Float64
          ωxc  = ωxp::Float64
          lac[2] += 1
        end

      #update ω1
      elseif up == 3

        ω1p = addupt(ω1c, ptn[3])::Float64

        # proposal likelihood and prior
        llr = 0.0

        # prior ratio
        prr = llrdbeta_x(ω1p, ω1c, ω1prior[1], ω1prior[2])

        if -Random.randexp()()() < (llr + prr)
          llc += llr::Float64
          prc += prr::Float64
          ω1c  = ω1p::Float64
          lac[3] += 1
        end

      # update ω0
      else

        ω0p = mulupt(ω0c, ptn[4])

        llr = 0.0

        # prior ratio
        prr = llrdexp_x(ω0p, ω0c, ω0prior)

        if -Random.randexp()()() < (llr + prr + log(ω0p/ω0c))
          llc += llr
          prc += prr
          ω0c  = ω0p
          lac[4] += 1
        end
      end

      # log number of updates
      ltn[up] += 1
      lup[up] += 1

      if (in(tune_int,ltn))
        wts = find(map(x -> x == tune_int,ltn))      # which to scale
        for j = wts
          ar     = lac[j]/lup[j]
          ptn[j] = scalef(ptn[j],ar)
          ltn[j] = 0
        end
      end
    end

    next!(p)
  end

  return llc, prc, σ²c, ωxc, ω1c, ω0c, λ1c, λ0c, ptn
end





