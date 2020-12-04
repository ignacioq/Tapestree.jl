#=

Anagenetic GBM birth-death MCMC

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#

function insane_gbmbd(tree    ::sTbd, 
                      out_file::String;
                      λprior  ::Float64           = 0.1,
                      μprior  ::Float64           = 0.1,
                      niter   ::Int64             = 1_000,
                      nthin   ::Int64             = 10,
                      nburn   ::Int64             = 200,
                      tune_int::Int64             = 100,
                      ϵi      ::Float64           = 0.4,
                      λi      ::Float64           = NaN,
                      μi      ::Float64           = NaN,
                      λtni    ::Float64           = 1.0,
                      μtni    ::Float64           = 1.0,
                      obj_ar  ::Float64           = 0.4,
                      pupdp   ::NTuple{4,Float64} = (0.4,0.4,0.1,0.1),
                      prints  ::Int64              = 5)

  # fix tree
  fixtree!(tree)

  # define δt
  th   = treeheight(tree)
  δt  *= th
  srδt = sqrt(δt)
  n    = sntn(tree)

   # starting parameters (using method of moments)
  if isnan(λi) && isnan(μi)
    λc, μc = moments(Float64(n), th, ϵi)
  else
    λc, μc = λi, μi
  end

  # make Ψ current and proposal parameters
  Ψc = iTgbmbd(tree, δt, srδt, log(λc), log(μc), σλi, σμi)
  Ψp = deepcopy(Ψc)

  # make fix Ψ directory
  idv = iDir[]
  bit = BitArray{1}()
  makeiDir!(Ψc, idv, bit)

  # make parent node directory to `iDir`
  inodes, terminus = make_inodes(idv)

  # create wbr vector
  wbr = falses(lastindex(idv))

  # create da branches vector
  dabr = Int64[]

  # make survival conditioning function (stem or crown)
  svf = iszero(pe(tree)) ? crown_prob_surv_cbd :
                           stem_prob_surv_cbd







