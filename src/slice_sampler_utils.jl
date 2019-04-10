#=

Slice sampling utilities

Ignacio Quintero Mächler

t(-_-t)

September 23 2017

=#





"""
    loop_slice_sampler(lhf      ::Function, 
                       p        ::Array{Float64,1},
                       nnps     ::Array{Int64,1},
                       nps      ::Array{Int64,1},
                       w        ::Array{Float64,1},
                       npars    ::Int64,
                       niter    ::Int64,
                       nthin    ::Int64)

Run slice sampling.
"""
function loop_slice_sampler(lhf      ::Function, 
                            p        ::Array{Float64,1},
                            nnps     ::Array{Int64,1},
                            nps      ::Array{Int64,1},
                            w        ::Array{Float64,1},
                            npars    ::Int64,
                            niter    ::Int64,
                            nthin    ::Int64)

  nlogs = fld(niter,nthin)
  its   = zeros(Float64,nlogs)
  hlog  = zeros(Float64,nlogs)
  ps    = zeros(Float64,nlogs,npars)

  lthin, lit = 0, 0

  # preallocate pp
  pp = copy(p)

  # start iterations
  prog = Progress(niter, 5, "running slice-sampler...", 20)

  hc = lhf(p)

  for it in Base.OneTo(niter) 

    for j in nnps
      S     = (hc - Random.randexp())
      L, R  = find_nonneg_int(p, pp, j, S, lhf, w[j])
      p, hc = sample_int(p, pp, j, L, R, S, lhf)
    end

    for j in nps
      S     = (hc - Random.randexp())
      L, R  = find_real_int(p, pp, j, S, lhf, w[j])
      p, hc = sample_int(p, pp, j, L, R, S, lhf)
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
    w_sampler(lhf, p::Array{Float64,1}, pupd::Array{Int64,1})

Run 100 iterations of the sampler to estimate appropriate w's.
"""
function w_sampler(lhf      ::Function, 
                   p        ::Array{Float64,1},
                   nnps     ::Array{Int64,1},
                   nps      ::Array{Int64,1},
                   npars    ::Int64,
                   optimal_w::Float64)

  w  = ones(npars)
  ps = Array{Float64,2}(undef, 100, npars)

  # posterior
  hc = lhf(p)

  # preallocate pp
  pp = Array{Float64,1}(undef, npars)

  prog = Progress(100, 5, "estimating optimal widths...", 20)

  for it in Base.OneTo(100)

    for j in nnps
     S     = (hc - Random.randexp())
     L, R  = find_nonneg_int(p, pp, j, S, lhf, w[j])
     p, hc = sample_int(p, pp, j, L, R, S, lhf)
    end

    for j in nps
      S     = (hc - Random.randexp())
      L, R  = find_real_int(p, pp, j, S, lhf, w[j])
      p, hc = sample_int(p, pp, j, L, R, S, lhf)
    end

    @inbounds setindex!(ps, p, it, :)

    next!(prog)
  end

  w = optimal_w .* (reduce(max, ps, dims=1) .- reduce(min, ps, dims=1))
  w = reshape(w, size(w,2))

  return (p, w)::Tuple{Array{Float64,1},Array{Float64,1}}
end





"""
    find_nonneg_int(p    ::Array{Float64}, 
                    pp   ::Array{Float64},
                    j    ::Int64, 
                    S    ::Float64, 
                    postf::Function, 
                    w    ::Float64)

Estimate a non_negative slice interval.
"""
function find_nonneg_int(p    ::Array{Float64}, 
                         pp   ::Array{Float64},
                         j    ::Int64, 
                         S    ::Float64, 
                         postf::Function, 
                         w    ::Float64)

  copyto!(pp, p)

  L::Float64 = pp[j] - w*rand()
  R::Float64 = L + w

  if L <= 0.0
    L = 1e-30
  end

  # left extreme
  pp[j] = L::Float64
  while S < postf(pp)
    L -= w::Float64
    if L <= 0.0
      L = 1e-30
      break
    end
    pp[j] = L::Float64
  end

  # right extreme
  pp[j] = R::Float64
  while S < postf(pp)
    R    += w::Float64
    pp[j] = R::Float64
  end

  return (L, R)::NTuple{2,Float64}
end





"""
    find_real_int(p    ::Array{Float64}, 
                  pp   ::Array{Float64}, 
                  j    ::Int64, 
                  S    ::Float64, 
                  postf::Function, 
                  w    ::Float64)

Estimate a non_negative slice interval.
"""
function find_real_int(p    ::Array{Float64}, 
                       pp   ::Array{Float64}, 
                       j    ::Int64, 
                       S    ::Float64, 
                       postf::Function, 
                       w    ::Float64)

  copyto!(pp, p)

  L::Float64 = pp[j] - w*rand()
  R::Float64 = L + w

  # left extreme
  pp[j] = L::Float64
  while S < postf(pp)
    L    -= w::Float64
    pp[j] = L::Float64
  end

  # right extreme
  pp[j] = R::Float64
  while S < postf(pp)
    R    += w::Float64
    pp[j] = R::Float64
  end

  return (L, R)::NTuple{2,Float64}
end





"""
    sample_int(p    ::Array{Float64,1}, 
               pp   ::Array{Float64,1},
               j    ::Int64, 
               L    ::Float64, 
               R    ::Float64, 
               S    ::Float64, 
               postf::Function)

Take one sample within the interval of the slice.
"""
function sample_int(p    ::Array{Float64,1}, 
                    pp   ::Array{Float64,1},
                    j    ::Int64, 
                    L    ::Float64, 
                    R    ::Float64, 
                    S    ::Float64, 
                    postf::Function)

  @inbounds begin
    copyto!(pp, p)


    while true
      pp[j] = (L + rand()*(R-L))::Float64

      postc = postf(pp)::Float64
      if S < postc
        copyto!(p, pp)
        return (p, postc)::Tuple{Array{Float64,1}, Float64}
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
    logdnorm(x::Float64, μ::Float64, σ²::Float64)
  
Compute the logarithmic transformation of the 
**Normal** density with mean `μ` and variance `σ²` for `x`.
"""
logdnorm(x::Float64, μ::Float64, σ²::Float64) = 
  @fastmath -(0.5*log(2.0π) + 0.5*log(σ²) + abs2(x - μ)/(2.0 * σ²))::Float64





"""
    logdexp(x::Float64, λ::Float64)

Compute the logarithmic transformation of the 
**Exponential** density with mean `λ` for `x`.
"""
logdexp(x::Float64, λ::Float64) = @fastmath (log(λ) - λ * x)::Float64





"""
    make_logdunif(a::Float64, b::Float64)

Make function to compute the logarithmic transformation of the 
**Uniform** density with lower bound `a` and upper bound `b` for `x`.
"""
function logdunif(x::Float64, a::Float64, b::Float64)
  @fastmath begin
    if x < a 
      return -Inf
    elseif x <= b 
      return -log(b-a)::Float64
    else 
      return -Inf
    end
  end
end





