#=

SSE richness distribution

Ignacio Quintero Mächler

t(-_-t)

23 07 2021
=#



@inline function pdt!(p1  ::Array{Float64,3}, 
                      p   ::Array{Float64,3}, 
                      λk  ::Float64, 
                      λj  ::Float64, 
                      λw  ::Float64, 
                      μk  ::Float64, 
                      μj  ::Float64,
                      gkj ::Float64,
                      gjk ::Float64,
                      nspp::Int64,
                      dt  ::Float64)

  @inbounds @fastmath begin

    for w in Base.OneTo(nspp - 1) 
      nw = Float64(w - 1)
      wh1 = w > 1

      for j in Base.OneTo(nspp - 1) 
        nj  = Float64(j - 1)
        jh1 = j > 1

        @simd for k in Base.OneTo(nspp - 1)
          nk = Float64(k - 1)
          kh1 = k > 1

          # if 0, 0, 0
          if !kh1 && !jh1 && !wh1
 
            p1[1,1,1] = 0.0

          else

            p1[k,j,w] += 
              (
                ## events
                # speciation in k
                (nk - 1.0 + nw) * λk * (kh1 ? p[k-1,j,w] : 0.0)         + 
                # speciation in j
                (nj - 1.0 + nw) * λj * (jh1 ? p[k,j-1,w] : 0.0)         + 
                # widespread speciation
                (nw + 1.0) * λw * ((jh1 && kh1) ? p[k-1,j-1,w+1] : 0.0) +
                # global extinction in k
                (nk + 1.0) * μk * p[k+1,j,w]                            + 
                # global extinction in j
                (nj + 1.0) * μj * p[k,j+1,w]                            +
                # local extinction in k
                (nw + 1.0) * μk * (jh1 ? p[k,j-1,w+1] : 0.0)            + 
                # local extinction in j
                (nw + 1.0) * μj * (kh1 ? p[k-1,j,w+1] : 0.0)            + 
                # colonization k -> j
                (nk + 1.0) * gkj * (wh1 ? p[k+1,j,w-1] : 0.0)           + 
                # colonization j -> k
                (nj + 1.0) * gjk * (wh1 ? p[k,j+1,w-1] : 0.0)           - 

                ## non events
                (
                  # speciation in k
                  (nk + nw) * λk + 
                  # speciation in j
                  (nj + nw) * λj +
                  # widespread speciation
                  nw * λw        +
                  # global extinction in k
                  nk * μk        + 
                  # global extinction in j
                  nj * μj        +
                  # local extinction in k
                  nw * μk        +
                  # local extinction in j
                  nw * μj        + 
                  # colonization k -> j
                  nk * gkj       + 
                  # colonization j -> k
                  nj * gjk
                ) * p[k,j,w]
              ) * dt
          end
        end
      end
    end
  end

  return nothing
end




"""
    _pdt!(t   ::Float64,
          p1  ::Array{Float64,3}, 
          p   ::Array{Float64,3}, 
          λk  ::Float64, 
          λj  ::Float64, 
          λw  ::Float64, 
          μk  ::Float64, 
          μj  ::Float64,
          gkj ::Float64,
          gjk ::Float64,
          nssp::Int64,
          dt  ::Float64)

Richness probabilities after time t.
"""
function _pdt!(t   ::Float64,
               p1  ::Array{Float64,3}, 
               p   ::Array{Float64,3}, 
               λk  ::Float64, 
               λj  ::Float64, 
               λw  ::Float64, 
               μk  ::Float64, 
               μj  ::Float64,
               gkj ::Float64,
               gjk ::Float64,
               nspp::Int64,
               dt  ::Float64)

  n = fld(t, dt)
  tl = nspp^3

  for i in Base.OneTo(Int64(n))
    pdt!(p1, p, λk, λj, λw, μk, μj, gkj, gjk, nspp, dt)
    unsafe_copyto!(p, 1, p1, 1, tl)
  end

  # last time
  m  = mod(t, dt)
  pdt!(p1, p, λk, λj, λw, μk, μj, gkj, gjk, nspp, m)

  return nothing
end




"""
    p_t(t   ::Float64,
        nspp::Int64,
        init::NTuple{3,Int64};
        dt  ::Float64 = 1e-3,
        λk  ::Float64 = 0.1,
        λj  ::Float64 = 0.1,
        λw  ::Float64 = 0.1,
        μk  ::Float64 = 0.01,
        μj  ::Float64 = 0.01,
        gkj ::Float64 = 0.1,
        gjk ::Float64 = 0.1)

Richness probabilities after time t given initial values and
number of species to consider.
"""
function p_t(t   ::Float64,
             nspp::Int64,
             init::NTuple{3,Int64};
             dt  ::Float64 = 1e-3,
             λk  ::Float64 = 0.1,
             λj  ::Float64 = 0.1,
             λw  ::Float64 = 0.1,
             μk  ::Float64 = 0.01,
             μj  ::Float64 = 0.01,
             gkj ::Float64 = 0.1,
             gjk ::Float64 = 0.1)

  p = zeros(nspp, nspp, nspp)
  p[init[1]+1, init[2]+1, init[3]+1] = 1.0
  p1 = copy(p)

  _pdt!(t, p1, p, λk, λj, λw, μk, μj, gkj, gjk, nspp, dt)

  return p1
end





