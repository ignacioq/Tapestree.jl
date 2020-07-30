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
                       p           ::Array{Array{Float64,1},1},
                       fp          ::Array{Array{Float64,1},1},
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
                       T           ::Float64,
                       screen_print::Int64)

  # estimate optimal w
  p, fp, w, o, t = 
    w_sampler(lhf, p, fp, nnps, nps, phid, mvps, nngps, mvhfs, hfgps,
      npars, optimal_w, screen_print, nburn, ntakew, nswap, ncch, winit, dt)



  # slice-sampler
  its, hlog, ps, Olog = 
    loop_slice_sampler(lhf, p, fp, nnps, nps, phid, mvps, nngps, mvhfs, hfgps, 
      w, npars, niter, nthin, nswap, ncch, o, t, screen_print)


  # choose cold chain (is equal to 1)
  P = Array{Float64,2}(undef, length(its), npars)
  H = Array{Float64,1}(undef, length(its))
  @inbounds begin
    for i in Base.OneTo(length(its))
      @views ii = findfirst(x -> isone(x), Olog[i,:])
      P[i,:] = ps[ii][i,:]
      H[i]   = hlog[i,ii]
    end
  end

  # save samples
  R = hcat(its, H, P)

  return R
end










