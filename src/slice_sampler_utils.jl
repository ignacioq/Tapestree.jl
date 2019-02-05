#=

Slice sampling utilities

Ignacio Quintero Mächler

t(-_-t)

September 23 2017

=#



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

  w = fill(1.0, npars)

  ps = zeros(Float64, 100, npars)

  hc::Float64 = lhf(p)

  prog = Progress(100, 5, "estimating optimal widths...", 20)

  for it in Base.OneTo(100)

    for j in nnps
      S     = (hc - Random.randexp())::Float64
      L, R  = find_nonneg_int(p, j, S, lhf, w[j])::NTuple{2,Float64}
      p, hc = sample_int(p, j, L, R, S, lhf)::Tuple{Array{Float64,1},Float64}
    end

    for j in nps
      S     = (hc - Random.randexp())::Float64
      L, R  = find_real_int(p, j, S, lhf, w[j])::NTuple{2,Float64}
      p, hc = sample_int(p, j, L, R, S, lhf)::Tuple{Array{Float64,1},Float64}
    end

    @inbounds setindex!(ps, p, it, :)

    next!(prog)
  end

  w = optimal_w .* (reduce(max, ps, dims=1) .- reduce(min, ps, dims=1))
  w = reshape(w, size(w,2))

  return (p, w)::Tuple{Array{Float64,1},Array{Float64,1}}
end






"""
    w_sampler(lhf, p::Array{Float64,1}, pupd::Array{Int64,1})

Run 100 iterations of the sampler to estimate appropriate w's.
"""
function w_sampler(lhf      ::Function, 
                   p        ::Array{Float64,1},
                   nnps     ::Array{Int64,1},
                   npars    ::Int64,
                   optimal_w::Float64)

  w = fill(1.0, npars)

  ps = zeros(Float64,100, npars)

  hc::Float64 = lhf(p)

  prog = Progress(100, 5, "estimating optimal w widths...", 20)

  for it in Base.OneTo(100)

    for j in nnps
      S     = (hc - Random.randexp())::Float64
      L, R  = find_nonneg_int(p, j, S, lhf, w[j])::NTuple{2,Float64}
      p, hc = sample_int(p, j, L, R, S, lhf)::Tuple{Array{Float64,1},Float64}
    end

    @inbounds setindex!(ps, p, it, :)

    next!(prog)
  end

  w = optimal_w .* (reduce(max, ps, dims=1) .- reduce(min, ps, dims=1))
  w = reshape(w, size(w,2))

  return (p, w)::Tuple{Array{Float64,1},Array{Float64,1}}
end







"""
    find_nonneg_int(p::Array{Float64}, j::Int64, S::Float64, postf, w::Float64)

Estimate a non_negative slice interval.
"""
function find_nonneg_int(p::Array{Float64}, j::Int64, S::Float64, postf, w::Float64)

  pc::Array{Float64,1} = deepcopy(p)

  L::Float64 = pc[j] - w*rand()
  R::Float64 = L + w

  if L <= 0.0
    L = 1e-30
  end

  # left extreme
  pc[j] = L::Float64
  while S < postf(pc)
    L -= w::Float64
    if L <= 0.0
      L = 1e-30
      break
    end
    pc[j] = L::Float64
  end

  # right extreme
  pc[j] = R::Float64
  while S < postf(pc)
    R    += w::Float64
    pc[j] = R::Float64
  end

  return (L, R)::NTuple{2,Float64}
end





"""
    find_real_int(p::Array{Float64}, j::Int64, S::Float64, postf, w::Float64)

Estimate a non_negative slice interval.
"""
function find_real_int(p::Array{Float64}, j::Int64, S::Float64, postf, w::Float64)

  pc::Array{Float64,1} = deepcopy(p)

  L::Float64 = pc[j] - w*rand()
  R::Float64 = L + w

  # left extreme
  pc[j] = L::Float64
  while S < postf(pc)
    L    -= w::Float64
    pc[j] = L::Float64
  end

  # right extreme
  pc[j] = R::Float64
  while S < postf(pc)
    R    += w::Float64
    pc[j] = R::Float64
  end

  return (L, R)::NTuple{2,Float64}
end






"""
    sample_int(p::Array{Float64}, j::Int64, L::Float64, R::Float64, S::Float64, postf)

Take one sample within the interval of the slice.
"""
function sample_int(p::Array{Float64}, 
                    j::Int64, 
                    L::Float64, 
                    R::Float64, 
                    S::Float64, 
                    postf)

  @inbounds begin

    pc::Array{Float64,1} = deepcopy(p)

    pc[j] = (L + rand()*(R-L))::Float64
    
    while true
      pc[j] = (L + rand()*(R-L))::Float64

      postc = postf(pc)::Float64
      if S < postc
        return (pc, postc)::Tuple{Array{Float64,1}, Float64}
      end

      if pc[j] < p[j]
        L = pc[j]::Float64
      else
        R = pc[j]::Float64
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





