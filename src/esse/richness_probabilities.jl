#=

SSE richness distribution

Ignacio Quintero Mächler

t(-_-t)

23 07 2021
=#




"""
    pdt_k!(pkt1::Array{Float64,1}, 
           pkt ::Array{Float64,1}, 
           pwt ::Array{Float64,1}, 
           λk  ::Float64, 
           λw  ::Float64, 
           μk  ::Float64, 
           μj  ::Float64, 
           gkj ::Float64,
           dt  ::Float64,
           lp  ::Int64)

One step for single area `k` in sse_g
"""
@inline function pdt_k!(pkt1::Array{Float64,1}, 
                        pkt ::Array{Float64,1}, 
                        pwt ::Array{Float64,1}, 
                        λk  ::Float64, 
                        λw  ::Float64, 
                        μk  ::Float64, 
                        μj  ::Float64, 
                        gkj ::Float64,
                        dt  ::Float64,
                        lp  ::Int64)

  @inbounds begin
    bs = ws = be = 0.0
    @simd for i in Base.OneTo(lp-1)
      pwti   = pwt[i]
      nw     = Float64(i - 1)
      bs  += nw * λw * pwti
      ws  += nw * λk * pwti
      be  += nw * μj * pwti
    end

    @simd for i in Base.OneTo(lp-1)
      pkti   = pkt[i]
      pktim1 = i > 1 ? pkt[i-1] : 0.0
      pktip1 = pkt[i+1]
      nk     = Float64(i - 1)

      pkt1[i] +=
        (nk + 1.0) * (μk + gkj)               * pktip1 * dt + 
        ((nk - 1.0) * λk + bs + ws + be)      * pktim1 * dt - 
        (nk * (λk + μk + gkj) + bs + ws + be) * pkti   * dt
    end
  end

  return nothing
end




"""
    pdt_w!(pwt1::Array{Float64,1}, 
           pwt ::Array{Float64,1}, 
           pkt ::Array{Float64,1}, 
           pjt ::Array{Float64,1}, 
           λw  ::Float64, 
           μk  ::Float64, 
           μj  ::Float64, 
           gkj ::Float64,
           gjk ::Float64,
           dt  ::Float64,
           lp  ::Int64)

One step for widespread area `w` in sse_g
"""
@inline function pdt_w!(pwt1::Array{Float64,1}, 
                        pwt ::Array{Float64,1}, 
                        pkt ::Array{Float64,1}, 
                        pjt ::Array{Float64,1}, 
                        λw  ::Float64, 
                        μk  ::Float64, 
                        μj  ::Float64, 
                        gkj ::Float64,
                        gjk ::Float64,
                        dt  ::Float64,
                        lp  ::Int64)

  @inbounds begin

    sk = sj = 0.0
    @simd for i in Base.OneTo(lp-1)
      n    = Float64(i - 1)
      sk  += n * gkj * pkt[i]
      sj  += n * gjk * pjt[i]
    end

    @simd for i in Base.OneTo(lp-1)
      pwti   = pwt[i]
      pwtim1 = i > 1 ? pwt[i-1] : 0.0
      pwtip1 = pwt[i+1]
      nw = Float64(i - 1)

      pwt1[i] +=
        (nw + 1.0) * (μk + μj + λw)     * pwtip1 * dt + 
        (sk + sj)                       * pwtim1 * dt -
        (nw * (μk + μj + λw) + sk + sj) * pwti   * dt 
    end
  end

  return nothing
end




"""
    _p_t!(t ::Float64,
          pkt1::Array{Float64,1}, 
          pjt1::Array{Float64,1}, 
          pwt1::Array{Float64,1}, 
          pkt ::Array{Float64,1}, 
          pjt ::Array{Float64,1}, 
          pwt ::Array{Float64,1}, 
          λk::Float64,
          λj::Float64,
          λw::Float64,
          μk::Float64,
          μj::Float64,
          gkj::Float64,
          gjk::Float64,
          dt ::Float64, 
          lp ::Int64)

Richness probabilities after time t.
"""
function _p_t!(t ::Float64,
               pkt1::Array{Float64,1}, 
               pjt1::Array{Float64,1}, 
               pwt1::Array{Float64,1}, 
               pkt ::Array{Float64,1}, 
               pjt ::Array{Float64,1}, 
               pwt ::Array{Float64,1}, 
               λk::Float64,
               λj::Float64,
               λw::Float64,
               μk::Float64,
               μj::Float64,
               gkj::Float64,
               gjk::Float64,
               dt ::Float64, 
               lp ::Int64)

  n = fld(t, dt)

  unsafe_copyto!(pkt1, 1, pkt, 1, lp)
  unsafe_copyto!(pjt1, 1, pjt, 1, lp)
  unsafe_copyto!(pwt1, 1, pwt, 1, lp)

  for i in Base.OneTo(Int64(n))
    pdt_k!(pkt1, pkt, pwt, λk, λw, μk, μj, gkj, dt, lp)
    pdt_k!(pjt1, pjt, pwt, λj, λw, μj, μk, gjk, dt, lp)
    pdt_w!(pwt1, pwt, pkt, pjt, λw, μk, μj, gkj, gjk, dt, lp)

    unsafe_copyto!(pkt, 1, pkt1, 1, lp)
    unsafe_copyto!(pjt, 1, pjt1, 1, lp)
    unsafe_copyto!(pwt, 1, pwt1, 1, lp)
  end

  # last time
  m  = mod(t, dt)
  pdt_k!(pkt1, pkt, pwt, λk, λw, μk, μj, gkj, m, lp)
  pdt_k!(pjt1, pjt, pwt, λj, λw, μj, μk, gjk, m, lp)
  pdt_w!(pwt1, pwt, pkt, pjt, λw, μk, μj, gkj, gjk, m, lp)

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

  pkt  = zeros(nspp)
  pjt  = zeros(nspp)
  pwt  = zeros(nspp)
  pkt1 = zeros(nspp)
  pjt1 = zeros(nspp)
  pwt1 = zeros(nspp)

  # starting values
  pkt[init[1]+1] = 1.0
  pjt[init[2]+1] = 1.0
  pwt[init[3]+1] = 1.0

  _p_t!(t, pkt1, pjt1, pwt1, pkt, pjt, pwt, 
    λk, λj, λw, μk, μj, gkj, gjk, dt, nspp)

  return pkt1, pjt1, pwt1
end





