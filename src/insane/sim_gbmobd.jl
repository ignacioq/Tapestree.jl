#=

Occurrence birth-death diffusion simulation

Jérémy Andréoletti

v(^-^v)

Created 22 09 2023
=#






#=
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Sample conditional on number of species
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#




"""
    sim_gbmobd(n       ::Int64;
               λ0      ::Float64         = 1.0,
               μ0      ::Float64         = 0.2,
               α       ::Float64         = 0.0,
               σλ      ::Float64         = 0.1,
               σμ      ::Float64         = 0.1,
               ψ       ::Vector{Float64} = [0.1],
               ω       ::Vector{Float64} = [1.0],
               ψts     ::Vector{Float64} = Float64[],
               init    ::Symbol          = :stem,
               δt      ::Float64         = 1e-3,
               nstar   ::Int64           = 2*n,
               p       ::Float64         = 5.0,
               warnings::Bool            = true,
               maxt    ::Float64         = δt*1e6)

Simulate `iTfbd` with `n` tips according to geometric Brownian motions for birth 
and death rates, piecewise-constant fossil sampling rates `ψ` (fossils included 
in the tree) and `ω`(fossil occurrences), with time shifts `ψts` between epochs.
"""
function sim_gbmobd(n       ::Int64;
                    λ0      ::Float64         = 1.0,
                    μ0      ::Float64         = 0.2,
                    α       ::Float64         = 0.0,
                    σλ      ::Float64         = 0.1,
                    σμ      ::Float64         = 0.1,
                    ψ       ::Vector{Float64} = [0.1],
                    ω       ::Vector{Float64} = [1.0],
                    ψts     ::Vector{Float64} = Float64[],
                    init    ::Symbol          = :stem,
                    δt      ::Float64         = 1e-3,
                    nstar   ::Int64           = 2*n,
                    p       ::Float64         = 5.0,
                    warnings::Bool            = true,
                    maxt    ::Float64         = δt*1e6)

  @assert issorted(ψts) "ψts should be sorted in increasing order"

  tree   = sim_gbmfbd(n, λ0=λ0, μ0=μ0, α=α, σλ=σλ, σμ=σμ, ψ=ψ, ψts=ψts, init=init, δt=δt, nstar=nstar, p=p, warnings=warnings, maxt=maxt)
  ωtimes = sim_occurrences(tree, ω, treeheight(tree) .- ψts)

  return tree, ωtimes
end


"""
    sim_occurrences(tree::T, ω::Vector{Float64}, ωts::Vector{Float64}) where {T <: iTree}

Simulate fossil occurrences on a given tree, with piecewise-constant rate `ω` and time shifts `ωts`.
"""
function sim_occurrences(tree::T, ω::Vector{Float64}, ωts::Vector{Float64}) where {T <: iTree}
    
    ωtimes = Float64[]
    e(tree)>1e-12 || def1(tree) || def2(tree) || return ωtimes
    LTT = ltt(tree)

    iω = 1
    ωc = ω[iω]

    while iω <= lastindex(ωts) && ωts[iω] >= LTT.t[1]
      iω += 1
      ωc = ω[iω]
    end

    for i in 1:(lastindex(LTT.t) - 1)
        t_start, t_end = LTT.t[i], LTT.t[i + 1]
        if !isapprox(t_start, t_end, atol=1e-12)
          ni = LTT.n[i]

          # Check if we have reached a fossilization rate shift
          while iω <= lastindex(ωts) && ωts[iω] >= t_end
              # Calculate the part of the interval before the rate shift
              Δt = t_start - ωts[iω]

              # Simulate occurrences using the Poisson process
              # @show t_start, ωts[iω]
              # @show ωc, ni, Δt
              nω = rand(Poisson(ωc * ni * Δt))
              ωtimes = append!(ωtimes, rand(Uniform(ωts[iω], t_start), nω))

              # Update values for the next segment
              t_start = ωts[iω]
              iω += 1
              ωc = ω[iω]
          end

          # Handle the remaining segment
          Δt = t_start - t_end
          # @show t_start, t_end
          # @show ωc, ni, Δt
          nω = rand(Poisson(ωc * ni * Δt))
          ωtimes = append!(ωtimes, rand(Uniform(t_end, t_start), nω))
        end
    end

    return ωtimes
end




# #=
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Sample conditional on time
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# =#





"""
    sim_gbmobd(t   ::Float64;
               λ0  ::Float64         = 1.0,
               μ0  ::Float64         = 0.2,
               α   ::Float64         = 0.0,
               σλ  ::Float64         = 0.1,
               σμ  ::Float64         = 0.1,
               ψ   ::Vector{Float64} = [0.1],
               ω   ::Vector{Float64} = [1.0],
               ψts ::Vector{Float64} = Float64[],
               δt  ::Float64         = 1e-3,
               nlim::Int64           = 10_000,
               init::Symbol          = :crown)

Simulate `iTfbd` of age `t` according to geometric Brownian motions for birth 
and death rates, piecewise-constant fossil sampling rates `ψ` (fossils included 
in the tree) and `ω`(fossil occurrences), with time shifts `ψts` between epochs.
"""
function sim_gbmobd(t   ::Float64;
                    λ0  ::Float64         = 1.0,
                    μ0  ::Float64         = 0.2,
                    α   ::Float64         = 0.0,
                    σλ  ::Float64         = 0.1,
                    σμ  ::Float64         = 0.1,
                    ψ   ::Vector{Float64} = [0.1],
                    ω   ::Vector{Float64} = [1.0],
                    ψts ::Vector{Float64} = Float64[],
                    δt  ::Float64         = 1e-3,
                    nlim::Int64           = 10_000,
                    init::Symbol          = :crown)

  @assert issorted(ψts, rev=true) "ψts should be sorted in decreasing order"

  tree   = sim_gbmfbd(t, λ0=λ0, μ0=μ0, α=α, σλ=σλ, σμ=σμ, ψ=ψ, ψts=ψts, δt=δt, nlim=nlim, init=init)
  ωtimes = sim_occurrences(tree, ω, ψts)

  return tree, ωtimes
end



