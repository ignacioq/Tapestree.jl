#=

GBM birth-death likelihood

Ignacio Quintero Mächler

t(-_-t)

Created 03 09 2020
=#




"""
    llik_gbm(tree::iTgbmbd, 
             σλ  ::Float64, 
             σμ  ::Float64,
             δt  ::Float64,
             srδt::Float64)

Returns the log-likelihood for a `iTgbmbd` according to GBM birth-death.
"""
function llik_gbm(tree::iTgbmbd, 
                  σλ  ::Float64, 
                  σμ  ::Float64,
                  δt  ::Float64,
                  srδt::Float64)

  tsb = ts(tree)
  lλb = lλ(tree)
  lμb = lμ(tree)

  if istip(tree) 
    ll_gbm_b(tsb, lλb, lμb, σλ, σμ, δt, srδt) + 
    (isextinct(tree) ? lμb[end] : 0.0)
  else
    ll_gbm_b(tsb, lλb, lμb, σλ, σμ, δt, srδt)    + 
    log(2.0) + lλb[end]                          + 
    llik_gbm(tree.d1::iTgbmbd, σλ, σμ, δt, srδt) + 
    llik_gbm(tree.d2::iTgbmbd, σλ, σμ, δt, srδt)
  end
end




"""
    ll_gbm_b(t   ::Array{Float64,1},
             lλv ::Array{Float64,1},
             lμv ::Array{Float64,1},
             σλ  ::Float64,
             σμ  ::Float64,
             δt  ::Float64, 
             srδt::Float64)

Returns the log-likelihood for a branch according to GBM birth-death.
"""
function ll_gbm_b(t   ::Array{Float64,1},
                  lλv ::Array{Float64,1},
                  lμv ::Array{Float64,1},
                  σλ  ::Float64,
                  σμ  ::Float64,
                  δt  ::Float64, 
                  srδt::Float64)

  @inbounds @fastmath begin

    # estimate standard `δt` likelihood
    nI = lastindex(t)-2

    llλ  = 0.0
    llμ  = 0.0
    llbd = 0.0
    lλvi = lλv[1]
    lμvi = lμv[1]
    @simd for i in Base.OneTo(nI)
      lλvi1 = lλv[i+1]
      lμvi1 = lμv[i+1]
      llλ  += (lλvi1 - lλvi)^2
      llμ  += (lμvi1 - lμvi)^2
      llbd += exp(0.5*(lλvi + lλvi1)) + exp(0.5*(lμvi + lμvi1)) 
      lλvi  = lλvi1
      lμvi  = lμvi1
    end

    # add to global likelihood
    ll = llλ*(-0.5/((σλ*srδt)^2)) - Float64(nI)*(log(σλ*srδt) + 0.5*log(2.0π)) + 
         llμ*(-0.5/((σμ*srδt)^2)) - Float64(nI)*(log(σμ*srδt) + 0.5*log(2.0π))
    # add to global likelihood
    ll -= llbd*δt

    # add final non-standard `δt`
    δtf   = t[nI+2] - t[nI+1]
    srδtf = sqrt(δtf)
    lλvi1 = lλv[nI+2]
    lμvi1 = lμv[nI+2]
    ll += ldnorm_bm(lλvi1, lλvi, srδtf*σλ)
    ll += ldnorm_bm(lμvi1, lμvi, srδtf*σμ)
    ll -= δtf*(exp(0.5*(lλvi + lλvi1)) + exp(0.5*(lμvi + lμvi1)))
  end

  return ll
end





"""
    pnorm(x::Float64, y::Float64, μ::Float64, σ::Float64)

Cumulative probability between `x` and `y`, with `x < y` for a 
**Normal** Distribution with mean `μ` and standard deviation `σ`.
"""
function pnorm(x::Float64, y::Float64, μ::Float64, σ::Float64)
  iσ = 1.0/(σ*sqrt(2.0))
  0.5*(erf((x - μ)*iσ, (y - μ)*iσ))
end





"""
estimate survival probability
"""

λ  = 0.1
μ  = 0.01
σλ = 0.2
σμ = 0.1
μ  = 0.01
δt = 1e-2
srδt = sqrt(δt)


et = eδt(et, λ, μ,δt)

function eδt(et::Float64, λ::Float64, μ::Float64, δt::Float64)
  μ*δt + 
  (1.0 - μ*δt)*(λ*δt)*et^2 + 
  (1.0 - λ*δt)*(1.0 - μ*δt)*et 
end






# divide the normal in 11 quantiles
qnorm(p::Float64) = sqrt(2.0) * erfinv(2.0*p - 1.0)

pnorm(x::Float64, μ::Float64, σ::Float64) = 
  0.5*(1.0 + erf((x - μ)/(σ*sqrt(2.0))))

pnorm(0.0, 0.0, 0.1) - pnorm(-0.1, 0.0, 0.1)

pnorm(0.1, 0.0, 0.1) - pnorm(0.0, 0.0, 0.1)


ldnorm_bm(0.1, 0.0, 0.2*δt)


λx = -10.0:0.05:10.0
μx = -10.0:0.05:10.0



# make transition probabilities
λy = [dnorm_bm(x, 0.0, σλ*srδt) for x in λx]
μy = [dnorm_bm(x, 0.0, σμ*srδt) for x in μx]


# reduce to non-zeros only
λr = findfirst(!iszero, λy):findlast(!iszero, λy)
μr = findfirst(!iszero, μy):findlast(!iszero, μy)

# λy = λy[λr]
# λx = λx[λr]
# μy = μy[μr]
# μx = μx[μr]






et = eδt(et, λ, μ, δt)


pnorm(-0.1, 0.0, σλ*srδt) - pnorm(-0.15, 0.0, σλ*srδt)
pnorm(-0.05, 0.0, σλ*srδt) - pnorm(-0.1, 0.0, σλ*srδt)
pnorm(0.05, 0.0, σλ*srδt) - pnorm(0.0, 0.0, σλ*srδt)

pnorm(0.0, 0.05, 0.0, σλ*srδt)



λx  = -1.01:0.02:1.01
λt1 = -1.00:0.02:1.00
μx = -1.01:0.02:1.01

collect(-1.01:0.02:1.01)

# make pnorms for λs
λpr = [pnorm(λx[i], λx[i+1], 0.0, σλ*srδt) for i in 1:length(λx)-1]
μpr = [pnorm(μx[i], μx[i+1], 0.0, σμ*srδt) for i in 1:length(μx)-1]

pnorm(0.0, 0.0, σλ*srδt) - pnorm(-0.01, 0.0, σλ*srδt)

λx  = -1.01:0.02:1.01
λt1 = -1.00:0.02:1.00
μx = -1.01:0.02:1.01

collect(-1.01:0.02:1.01)

# make pnorms for λs
λpr = [pnorm(λx[i], λx[i+1], 0.0, σλ*srδt) for i in 1:length(λx)-1]
μpr = [pnorm(μx[i], μx[i+1], 0.0, σμ*srδt) for i in 1:length(μx)-1]

pnorm(0.0, 0.0, σλ*srδt) - pnorm(-0.01, 0.0, σλ*srδt)






"""
pnorm is symmetrical, only do half
"""


# 1D = lambda
# 2D = mu
# 3D = time
E[2,,] += 

nr, nc, nt = size(E)

ei = 1 # lambda dimension
ej = 1 # mu dimension
ek = 2 # time dimension

# make integrating function to sum over previous multiplied by probabilities
E[ei,ej,ek] =

ens = 0.0 
for i in 1:(nr-1), j in 1:(nc-1)
  ens += (λpr[i]*μpr[j]*
    "make the right indexes here"
    E[i,j,ek-1])
end

μ*δt + 
(1.0 - μ*δt)*(λ*δt) * (gλ[1]*E[1,,]*gμ[1]*E[1,,])^2 + 
(1.0 - λ*δt)*(1.0 - μ*δt) * (gλ[1]*E[1,,]*gμ[1]*E[1,,])


# indexes for probabilities
E1 = 51:101
E2 = 50:100
E3 = 49:99
E4 = 48:98
E5 = 47:97
E101 = 2:52
E102 = 1:51






"""
estimate survival probability
"""

λi  = 0.1
μi  = 0.01
σλ = 0.2
σμ = 0.1
δt = 1e-2
srδt = sqrt(δt)

tt  = 1.0
tis = 0.0:δt:tt

# this could be relative to the standard deviation
δN = 0.02 
Ndev = 1.0 

# make lambda values vector centered at λi 
lλ1 = collect((log(λi) - Ndev):δN:(log(λi) + Ndev))
λ1  = exp.(lλ1)
lμ1 = collect((log(μi) - Ndev):δN:(log(μi) + Ndev))
μ1  = exp.(lμ1)

# estimate probabilitites
λpr = [pnorm(i, i+δN, 0.0, σλ*srδt) 
  for i in (-2.0*Ndev - δN*0.5):δN:(2.0*Ndev - δN*0.5)]
μpr = [pnorm(i, i+δN, 0.0, σμ*srδt) 
  for i in (-2.0*Ndev - δN*0.5):δN:(2.0*Ndev - δN*0.5)]




λpr[101:201]
λpr[100:200]
.
.
.
λpr[1:101]

# make E matrix
# 1D = lambda
# 2D = mu
# 3D = time
E = zeros(length(λt1), length(μt1), length(tis))

nr, nc, nt = size(E)

ei = 1
ej = 1
et = 2


# function integrating over all 
for et in 2:nt, ej in Base.OneTo(nc), ei in Base.OneTo(nr)
  ens = 0.0 
  for j in Base.OneTo(nc), i in Base.OneTo(nr)
    ens += λpr[nr-ei+i] * μpr[nc-ej+j] * E[i,j,et-1]
  end

  E[ei,ej,et] = μ1[ej]*δt                         + 
    (1.0 - μ1[ej]*δt) * (λ1[ei]*δt)       * ens^2 + 
    (1.0 - λ1[ei]*δt) * (1.0 - μ1[ej]*δt) * ens
end


# now make function considering only the starting point
λr = findfirst(!iszero, λpr):findlast(!iszero, λpr)
μr = findfirst(!iszero, μpr):findlast(!iszero, μpr)


for et in 2:nt, ej in Base.OneTo(nc), ei in Base.OneTo(nr)
  ens = 0.0 
  for j in Base.OneTo(nc), i in Base.OneTo(nr)
    ens += λpr[nr-ei+i] * μpr[nc-ej+j] * E[i,j,et-1]
  end

  E[ei,ej,et] = E[ei,ej,et-1]                                 +
                μ1[ej]*δt                                     + 
                (1.0 - μ1[ej]*δt) * (λ1[ei]*δt)       * ens^2 + 
                (1.0 - λ1[ei]*δt) * (1.0 - μ1[ej]*δt) * ens
end



λr1 = λr[1]
λr2 = λr[end]



nr - 1 + 1 

"""
This could be estimated before hand just once
"""
ens = 0.0 
for j in Base.OneTo(nc), i in Base.OneTo(nr)
  ens += λpr[nr-ei+i] * μpr[nc-ej+j] * E[i,j,et-1]
end

E[51,51,2] = μ1[ej]*δt                         + 
    (1.0 - μ1[ej]*δt) * (λ1[ei]*δt)       * ens^2 + 
    (1.0 - λ1[ei]*δt) * (1.0 - μ1[ej]*δt) * ens



ens = 0.0 
for j in Base.OneTo(nc), i in Base.OneTo(nr)
  ens += λpr[nr-ei+i] * μpr[nc-ej+j] * E[i,j,et-1]
end

E[51,51,3] = 


λr
μr

(1:1,1:1)


"""
Simplify everything to have only three targets (left - middle - right).
"""








