#=

Slice sampling for EGeoHiSSE

Ignacio Quintero Mächler

t(-_-t)

November 20 2017

=#





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

Run slice-sampling Markov Chain given posterior function.
"""
function slice_sampler(lhf         ::Function, 
                       p           ::Array{Float64,1},
                       fp          ::Array{Float64,1},
                       nnps        ::Array{Int64,1},
                       nps         ::Array{Int64,1},
                       phid        ::Array{Int64,1},
                       mvps        ::Array{Array{Int64,1},1},
                       nngps       ::Array{Array{Bool,1},1},
                       mvhfs       ::Array{Array{Int64,1},1},
                       npars       ::Int64,
                       niter       ::Int64,
                       nthin       ::Int64,
                       nburn       ::Int64,
                       ntakew      ::Int64,
                       winit       ::Float64,
                       optimal_w   ::Float64,
                       screen_print::Int64)

  # estimate optimal w
  p, fp, w = 
    w_sampler(lhf, p, fp, nnps, nps, phid, mvps, nngps, mvhfs, 
      npars, optimal_w, screen_print, nburn, ntakew, winit)

  # slice-sampler
  its, hlog, ps = 
    loop_slice_sampler(lhf, p, fp, nnps, nps, phid, mvps, nngps, mvhfs, 
      w, npars, niter, nthin, screen_print)

  # save samples
  R = hcat(its, hlog, ps)

  return R
end





# """
#     slice_sampler(tip_val    ::Dict{Int64,Array{Float64,1}},
#                   edges      ::Array{Int64,2},
#                   edlen      ::Array{Float64,1},
#                   out_file   ::String;
#                   constraints::String  = NaN,
#                   niter      ::Int64   = 10_000,
#                   nthin      ::Int64   = 10,
#                   model      ::String  = "musse",
#                   λpriors    ::Float64 = .1,
#                   μpriors    ::Float64 = .1,
#                   qpriors    ::Float64 = .1)

# Run slice-sampling Markov Chain for MuSSE model.
# """
# function slice_sampler(tip_val    ::Dict{Int64,Array{Float64,1}},
#                        edges      ::Array{Int64,2},
#                        edlen      ::Array{Float64,1},
#                        out_file   ::String;
#                        constraints::NTuple{N,String} = (" ",),
#                        niter      ::Int64            = 10_000,
#                        nthin      ::Int64            = 10,
#                        λpriors    ::Float64          = .1,
#                        μpriors    ::Float64          = .1,
#                        qpriors    ::Float64          = .1) where {N}

#   # k areas
#   k = length(tip_val[1])

#   # make ode and define parameters
#   mod_ode = make_musse(k)

#   npars = k + k*k 
#   p = fill(0.2,npars)
#   δ = length(tip_val)/sum(edlen)
#   p[1:k]      .= δ + rand()    # set λs
#   p[(k+1):2k] .= p[1] - δ      # set μs

#   pardic = build_par_names(k)

#   # parameter update
#   pupd = Base.OneTo(npars)

#   #constraints
#   conp = set_constraints(constraints, pardic)

#   # remove contraints for being updated
#   pupd = setdiff(pupd, keys(conp))

#   # create likelihood, prior and posterior functions
#   llf = make_llf(tip_val, edges, edlen, mod_ode, p, sbrlen = sum(edlen))
#   lpf = make_lpf(λpriors, μpriors, qpriors, k)
#   lhf = make_lhf(llf, lpf, conp)

#   # set up slice-sampling
#   nlogs = fld(niter,nthin)
#   its   = zeros(Float64,nlogs)
#   h     = zeros(Float64,nlogs)
#   ps    = zeros(Float64,nlogs,npars)

#   lthin = 0::Int64
#   lit   = 0::Int64

#   # cat
#   printstyled("running MuSSE model \n", color=:green)

#   # estimate w
#   p, w = w_sampler(lhf, p, pupd, npars)

#   # start iterations
#   prog = Progress(niter, 5, "running slice-sampler....", 20)

#   hc = lhf(p)::Float64

#   for it in Base.OneTo(niter) 

#     for j in pupd
#       S     = hc - Random.randexp()
#       L, R  = find_nonneg_int(p, j, S, lhf, w[j])
#       p, hc = sample_int(p, j, L, R, S, lhf)
#     end

#     # log samples
#     lthin += 1
#     if lthin == nthin
#       @inbounds begin
#         lit += 1
#         setindex!(its, it, lit)
#         setindex!(h,   hc, lit)
#         setindex!(ps,   p, lit, :)
#       end
#       lthin = 0
#     end

#     next!(prog)
#   end

#   # save samples
#   R = hcat(its, h, ps)

#   # column names
#   col_nam = ["Iteration", "Posterior"]

#   # add λ names
#   for i in 0:(k-1)
#     push!(col_nam, "lamdba$i")
#   end

#   # add μ names
#   for i in 0:(k-1)
#     push!(col_nam, "mu$i")
#   end

#   # add q names
#   for j in 0:(k-1), i in 0:(k-1)
#     if i == j 
#       continue
#     end
#     push!(col_nam, "q$i$j")
#   end

#   R = vcat(reshape(col_nam, 1, lastindex(col_nam)), R)

#   writedlm(out_file*".log", R)

#   return R
# end




