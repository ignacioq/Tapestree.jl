#=

Model definition for *SSE models

Ignacio Quintero Mächler

t(-_-t)

27 05 2019

=#





"""
    define_mod(esse_mod::String, x::Array{Float64,1}, y::Array{Float64,1}, k::Int64)

Defines ESSE model.
"""
function define_mod(esse_mod::String,
                    k       ::Int64,
                    af      ::Function,
                    md      ::Bool)

  if occursin(r"^[s|S][A-za-z]*", esse_mod)         # if speciation
    mod_ode = md ? make_esse_s(k, af, md) : make_esse_s(k, af)
    npars   = 2k + k*k
    pardic  = build_par_names(k, true)
    ws      = true
    printstyled("running speciation ESSE model \n", color=:green)

  elseif occursin(r"^[e|E][A-za-z]*", esse_mod)     # if extinction
    mod_ode = md ? make_esse_e(k, af, md) : make_esse_e(k, af)
    npars   = (2k + k*k)::Int64
    pardic  = build_par_names(k, true)
    ws      = false

    printstyled("running extinction ESSE model \n", color=:green)

  elseif occursin(r"^[t|T|r|R|q|Q][A-za-z]*", esse_mod) # if transition
    mod_ode = md ? make_esse_q(k, af, md) : make_esse_q(k, af)
    npars   = 2k*k::Int64
    pardic  = build_par_names(k, false)
    ws      = false

    printstyled("running transition ESSE model \n", color=:green)

  else 
    error("esse_mod does not match any of the alternatives: 
          speciation, extinction or transition")
  end

  return mod_ode, npars, pardic, md, ws

end





"""
    define_mod(cov_mod::String,
               k            ::Int64,
               h            ::Int64,
               ny           ::Int64)

Defines EGeoHiSSE model for `k` areas, `h` hidden states and `ny` covariates.
"""
function define_mod(cov_mod::NTuple{N,String},
                    k      ::Int64,
                    h      ::Int64,
                    ny     ::Int64) where {N}

  model = [false, false, false]

  for m in cov_mod
    # if speciation
    if occursin(r"^[s|S][A-za-z]*", m) 
      model[1] = true
    end
    # if extinction
    if occursin(r"^[e|E][A-za-z]*", m) 
      model[2] = true
    end
    # if dispersal
    if occursin(r"^[t|T|r|R|q|Q][A-za-z]*", m)
      model[3] = true
    end
  end

  mexp = "$(model[1] ? "speciation," : "")$(model[2] ? "extinction," : "")$(model[3] ? "transition," : "")"
  mexp = replace(mexp, "," => ", ")
  mexp = mexp[1:(end-2)]

  printstyled("running $mexp EGeoHiSSE model with:
    $k single areas 
    $h hidden states 
    $ny covariates \n", color=:green)

  return tuple(model...)
end



