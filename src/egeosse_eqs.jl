#=

EGEOSSE equations

Ignacio Quintero Mächler

t(-_-t)

Created 15 02 2019
=#






"""
 "lambda_A" => 1
 "lambda_B" => 2
 "lambda_C" => 3
 "lambda_W" => 4
     "mu_A" => 5
     "mu_B" => 6
     "mu_C" => 7
  "gain_AB" => 8
  "gain_AC" => 9
  "gain_BA" => 10
  "gain_BC" => 11
  "gain_CA" => 12
  "gain_CB" => 13
   "loss_A" => 14
   "loss_B" => 15
   "loss_C" => 16
   "beta_A" => 17
   "beta_B" => 18
   "beta_C" => 19

    D"A"   => 1
    D"B"   => 2
    D"C"   => 3
    D"AB"  => 4
    D"AC"  => 5
    D"BC"  => 6
    D"ABC" => 7
    E"A"   => 8
    E"B"   => 9
    E"C"   => 10
    E"AB"  => 11
    E"AC"  => 12
    E"BC"  => 13
    E"ABC" => 14
"""

"""
    make_geosse(k::Int64)

GeoSSE ODE equation for k areas.
"""
function make_geosse(k::Int64)

  # n states
  ns = 2^k - 1

  #create subsets
  sa = Array{String,1}(undef, k)
  for i in Base.OneTo(k)
    sa[i] = string('A' + (i-1))
  end
  S = subsets(sa)
  popfirst!(S)            # remove empty
  sort!(S, by = length)   # arrange by range size

  # Metaprogram to generate GeoSSE equations
  eqs = quote end

  for ri = Base.OneTo(ns)

    # range
    r = S[ri]

    # length of r
    lr = lastindex(r)

    # which single areas occur in r
    ia = findall(x -> occursin(x, r), sa)
    oa = setdiff(1:k, ia)

    #= 
    likelihoods
    =#

    # no events
    nev = noevents_expr(ri, lr, ia, oa, k)

    # local extinction
    # remove if !isone(lr)
    lex = localext_expr(r, ia, sa, S, k)

    # dispersal
    # remove if lr == k
    dis = dispersal_expr(r, lr, ia, oa, S, k)

    # within-region speciation
    wrs = wrspec_expr(ri, ia, ns)

    # between-region speciation
    # remove if !isone(lr)
    brs = brspec_expr(r, S, ns)

    # if single area
    if ri <= k
      push!(eqs.args, 
        quote du[$ri] = $nev + $dis + $wrs end)
    # if widespread
    elseif ri != ns
      push!(eqs.args, 
        quote du[$ri] = $nev + $lex + $dis + $wrs + $brs end)
    else
      push!(eqs.args, 
        quote du[$ri] = $nev + $lex + $wrs + $brs end)
    end

    #= 
    extinctions
    =#


  end
  

end





"""
    noevents_expr(ri::Int64,
                  lr::Int64,
                  ia::Array{Int64,1},
                  oa::Array{Int64,1},
                  k ::Int64)

Return expression for no events.
"""
function noevents_expr(ri::Int64,
                       lr::Int64,
                       ia::Array{Int64,1},
                       oa::Array{Int64,1},
                       k ::Int64)

  ts = isone(lr) ? 0 : (k*(k-1) + k)

  ex = :(+ ($(2^(lr-1) - 1.) * p[$(k+1)]))
  for (i, a) = enumerate(ia)
    push!(ex.args, :(p[$a] + p[$(a + k + ts + 1)]))
    for j = oa
      j -= a <= j ? 1 : 0
      push!(ex.args[i+2].args, :(p[$(2k + 1 + (k-1)*(a-1) + j)]))
    end
  end
  ex = :(-1.0 * $ex * u[$ri])

  # remove 0 product if single area
  if isone(lr)
    ex.args[3] = ex.args[3].args[3]
  end

  return ex
end





"""
    localext_expr(r ::String,
                  ia::Array{Int64,1},
                  sa::Array{String,1},
                  S ::Array{String,1},
                  k ::Int64)

Return expression for local extinction.
"""
function localext_expr(r ::String,
                       ia::Array{Int64,1},
                       sa::Array{String,1},
                       S ::Array{String,1},
                       k ::Int64)

  ex = :(1+1)
  for a = ia
    push!(ex.args, :(p[$(k^2 + k + 1 + a)] * 
                     u[$(findfirst(x -> x == replace(r, sa[a] => ""), S))]))
  end
  deleteat!(ex.args, 2:3)

  return ex
end





"""
    dispersal_expr(r ::String,
                   lr::Int64,
                   ia::Array{Int64,1},
                   oa::Array{Int64,1},
                   S ::Array{String,1},
                   k ::Int64)

Return expression for dispersal.
"""
function dispersal_expr(r ::String,
                        lr::Int64,
                        ia::Array{Int64,1},
                        oa::Array{Int64,1},
                        S ::Array{String,1},
                        k ::Int64)

  ida = findall(x -> all(occursin.(split(r,""),x)) && 
                     lastindex(x) == (lr + 1), 
                S)
  ex = :(1+1)
  for a = ia, (i, j) = enumerate(oa)
    j -= a <= j ? 1 : 0
    push!(ex.args, :(p[$(2k + 1 + (k-1)*(a-1) + j)] * u[$(ida[i])]))
  end

  deleteat!(ex.args, 2:3)

  return ex
end





"""
    wrspec_expr(ri::Int64,
                ia::Array{Int64,1},
                ns::Int64)

Return expression for within-region speciation.
"""
function wrspec_expr(ri::Int64,
                     ia::Array{Int64,1},
                     ns::Int64)

  if isone(lastindex(ia)) 
    wrs = :(2.0 * p[$ri] * u[$(ri + ns)] * u[$ri])
  else
    wrs = :(1 + 1)
    for i = ia
      push!(wrs.args, :(p[$i] * (u[$(i + ns)] * u[$ri] + u[$(ri + ns)] * u[$i])))
    end
    deleteat!(wrs.args, 2:3)
  end

  return wrs
end




"""
    brspec_expr(r::String,
                S ::Array{String,1},
                ns::Int64)

Return expression for within-region speciation.
"""
function brspec_expr(r ::String,
                     S ::Array{String,1},
                     ns::Int64)

  va  = vicsubsets(r)
  brs = :(1+1)
  for (la, ra) = va
    push!(brs.args,
      :(u[$(findfirst(isequal(ra), S) + ns)] *
        u[$(findfirst(isequal(la), S))]))
  end
  deleteat!(brs.args, 2:3)
  brs = :($(2^lastindex(r) - 3.0) * p[$(k+1)] * $brs)

  isone(brs.args[2]) && deleteat!(brs.args, 2)

  return brs
end




brspec_expr("AB", S, 7)



brspec_expr("AC", S, 7)
brspec_expr("BC", S, 7)

brspec_expr("ABC", S, 7)


"lambda_A" => 1
"lambda_B" => 2
"lambda_C" => 3
"lambda_W" => 4
    "mu_A" => 5
    "mu_B" => 6
    "mu_C" => 7
 "gain_AB" => 8
 "gain_AC" => 9
 "gain_BA" => 10
 "gain_BC" => 11
 "gain_CA" => 12
 "gain_CB" => 13
  "loss_A" => 14
  "loss_B" => 15
  "loss_C" => 16
  "beta_A" => 17
  "beta_B" => 18
  "beta_C" => 19


  D"A"   => 1
  D"B"   => 2
  D"C"   => 3
  D"AB"  => 4
  D"AC"  => 5
  D"BC"  => 6
  D"ABC" => 7
  E"A"   => 8
  E"B"   => 9
  E"C"   => 10
  E"AB"  => 11
  E"AC"  => 12
  E"BC"  => 13
  E"ABC" => 14











  function f(du::Array{Float64,1}, 
             u::Array{Float64,1}, 
             p::Array{Float64,1}, 
             t::Float64)

  for j = Base.OneTo(ns)

    # no events
    for j = Base.OneTo(ns)


    du[j]    = 


    du[j+ns] =
  end





  function f(du::Array{Float64,1}, 
             u::Array{Float64,1}, 
             p::Array{Float64,1}, 
             t::Float64)

    @inbounds begin
      # probabilities
      du[1] = -(p[1] + p[4] + p[6])*u[1] + p[6]*u[3] + 2.0*p[1]*u[1]*u[4]
      du[2] = -(p[2] + p[5] + p[7])*u[2] + p[7]*u[3] + 2.0*p[2]*u[2]*u[5]
      du[3] = -(p[1] + p[2] + p[3] + p[4] + p[5])*u[3] + 
                p[4]*u[2] + p[5]*u[1]        +
                p[1]*(u[4]*u[3] + u[6]*u[1]) + 
                p[2]*(u[5]*u[3] + u[6]*u[2]) + 
                p[3]*(u[4]*u[2] + u[5]*u[1])

      # extinction
      du[4] = -(p[1] + p[4] + p[6])*u[4] + p[4] + p[6]*u[6] + p[1]*u[4]^2
      du[5] = -(p[2] + p[5] + p[7])*u[5] + p[5] + p[7]*u[6] + p[2]*u[5]^2
      du[6] = -(p[1] + p[2] + p[3] + p[4] + p[5])*u[6] + 
                p[4]*u[5] + p[5]*u[4] + p[1]*u[6]*u[4] + 
                p[2]*u[6]*u[5] + p[3]*u[4]*u[5]
    end

    return nothing
  end

  return f
end






#=
  Model with separate endemic extinction rates.

      1    2    3     4    5    6         7         8         9       
  p = sa1, sb1, sab1, ma1, mb1, qa1->ab1, qb1->ab1, qab1->a1, qab1->b1,

      10   11   12    13   14   15        16        17        18      
      sa2, sb2, sab2, ma2, mb2, qa2->ab2, qb2->ab2, qab2->a2, qab2->b2,

      19       20       21   22   23   24
      wi1->i2, wi2->i1, βa1, βb1, βa2, βb2,

      1    2    3     4    5    6     7    8    9     10   11   12
  u = da1, db1, dab1, ea1, eb1, eab1, da2, db2, dab2, ea2, eb2, eab2 
=#
"""
    egeohisse_2k_s(du::Array{Float64,1}, 
                   u::Array{Float64,1}, 
                   p::Array{Float64,1}, 
                   t::Float64)

Speciation EGeoHiSSE equation for 2 areas and 2 hidden states.
"""
function make_egeohisse_2k_s(af::Function)

  function f(du::Array{Float64,1}, 
             u::Array{Float64,1}, 
             p::Array{Float64,1}, 
             t::Float64)

    # UNIdimensional time function
    aft = af(t)

    sa1 = p[1]  * exp(aft*p[21])
    sb1 = p[2]  * exp(aft*p[22])
    sa2 = p[10] * exp(aft*p[23])
    sb2 = p[11] * exp(aft*p[24])

    @inbounds begin

      ### Hidden State 1
      ## probabilities 
      # area A
      du[1] = -(sa1 + p[4] + p[6] + p[19])*u[1]            + 
                p[6]*u[3] + p[19]*u[7] + 2.0*sa1*u[1]*u[4]
      # area B
      du[2] = -(sb1 + p[5] + p[7] + p[19])*u[2]            + 
                p[7]*u[3] + p[19]*u[8] + 2.0*sb1*u[2]*u[5]
      # area AB
      du[3] = -(sa1 + sb1 + p[3] + p[8] + p[9] + p[19])*u[3] + 
                p[8]*u[1] + p[9]*u[2]+ p[19]*u[9]            +
                sa1*(u[4]*u[3] + u[6]*u[1])                  + 
                sb1*(u[5]*u[3] + u[6]*u[2])                  + 
                p[3]*(u[4]*u[2] + u[5]*u[1])
      ## extinction
      # area A
      du[4] = -(sa1 + p[4] + p[6] + p[19])*u[4]             + 
                p[4] + p[19]*u[10] + p[6]*u[6] + sa1*u[4]^2
      # area B
      du[5] = -(sb1 + p[5] + p[7] + p[19])*u[5]             + 
                p[5] + p[19]*u[11] + p[7]*u[6] + sb1*u[5]^2
      # area AB
      du[6] = -(sa1 + sb1 + p[3] + p[8] + p[9] + p[19])*u[6]   +
                p[8]*u[4] + p[9]*u[5] + p[19]*u[12]              + 
                sa1*u[6]*u[4] + sb1*u[6]*u[5] + p[3]*u[4]*u[5]
      ### Hidden State 2
      ## probabilities 
      # area A
      du[7] = -(sa2 + p[13] + p[15] + p[20])*u[7]            + 
                p[15]*u[9] + p[20]*u[1] + 2.0*sa2*u[7]*u[10]
      # area B
      du[8] = -(sb2 + p[14] + p[16] + p[20])*u[8]            + 
                p[16]*u[9] + p[20]*u[2] + 2.0*sb2*u[8]*u[11]
      # area AB
      du[9] = -(sa2 + sb2 + p[12] + p[17] + p[18] + p[20])*u[9] + 
                p[17]*u[7] + p[18]*u[8]+ p[20]*u[3]                 +
                sa2*(u[10]*u[9] + u[12]*u[7])                     + 
                sb2*(u[11]*u[9] + u[12]*u[8])                     + 
                p[12]*(u[10]*u[8] + u[11]*u[7])
      ## extinction
      # area A
      du[10] = -(sa2 + p[13] + p[15] + p[20])*u[10]            + 
                 p[13] + p[20]*u[4] + p[15]*u[12] + sa2*u[10]^2
      # area B
      du[11] = -(sb2 + p[14] + p[16] + p[20])*u[11]            + 
                 p[14] + p[20]*u[5] + p[16]*u[12] + sb2*u[11]^2
      # area AB
      du[12] = -(sa2 + sb2 + p[12] + p[17] + p[18] + p[20])*u[6]       +
                 p[17]*u[10] + p[18]*u[11] + p[20]*u[6]                    + 
                 sa2*u[10]*u[12] + sb2*u[11]*u[12] + p[12]*u[10]*u[11]
    end

    return nothing
  end

  return f
end






#=
    1   2   3    4   5   6   7
p = sa, sb, sab, xa, xb, qa, qb

    1   2   3    4   5   6
u = da, db, dab, ea, eb, eab
=#
"""
    geosse_2k(du::Array{Float64,1}, 
              u::Array{Float64,1}, 
              p::Array{Float64,1}, 
              t::Float64)

GeoSSE ODE equation for 2 areas.
"""
function geosse_2k(du::Array{Float64,1}, 
                   u::Array{Float64,1}, 
                   p::Array{Float64,1}, 
                   t::Float64)

  @inbounds begin
    # probabilities
    du[1] = -(p[1] + p[4] + p[6])*u[1] + p[6]*u[3] + 2.0*p[1]*u[1]*u[4]
    du[2] = -(p[2] + p[5] + p[7])*u[2] + p[7]*u[3] + 2.0*p[2]*u[2]*u[5]
    du[3] = -(p[1] + p[2] + p[3] + p[4] + p[5])*u[3] + 
              p[4]*u[2] + p[5]*u[1]        +
              p[1]*(u[4]*u[3] + u[6]*u[1]) + 
              p[2]*(u[5]*u[3] + u[6]*u[2]) + 
              p[3]*(u[4]*u[2] + u[5]*u[1])

    # extinction
    du[4] = -(p[1] + p[4] + p[6])*u[4] + p[4] + p[6]*u[6] + p[1]*u[4]^2
    du[5] = -(p[2] + p[5] + p[7])*u[5] + p[5] + p[7]*u[6] + p[2]*u[5]^2
    du[6] = -(p[1] + p[2] + p[3] + p[4] + p[5])*u[6] + 
              p[4]*u[5] + p[5]*u[4] + p[1]*u[6]*u[4] + 
              p[2]*u[6]*u[5] + p[3]*u[4]*u[5]
  end
  return nothing
end






#=
  Model with separate endemic extinction rates.

      1    2    3     4    5    6         7         8         9       
  p = sa1, sb1, sab1, ma1, mb1, qa1->ab1, qb1->ab1, qab1->a1, qab1->b1,

      10   11   12    13   14   15        16        17        18      
      sa2, sb2, sab2, ma2, mb2, qa2->ab2, qb2->ab2, qab2->a2, qab2->b2,

      19       20     
      wi1->i2, wi2->i1

      1    2    3     4    5    6     7    8    9     10   11   12
  u = da1, db1, dab1, ea1, eb1, eab1, da2, db2, dab2, ea2, eb2, eab2 
=#
"""
    geohisse_2k(du::Array{Float64,1}, 
                u::Array{Float64,1}, 
                p::Array{Float64,1}, 
                t::Float64)

GeoHiSSE + extinction most general ODE equation for 2 areas & 2 hidden states.
"""
function geohisse_2k(du::Array{Float64,1}, 
                     u::Array{Float64,1}, 
                     p::Array{Float64,1}, 
                     t::Float64)

  @inbounds begin

    ### Hidden State 1
    ## probabilities 
    # area A
    du[1] = -(p[1] + p[4] + p[6] + p[19])*u[1]            + 
              p[6]*u[3] + p[19]*u[7] + 2.0*p[1]*u[1]*u[4]
    # area B
    du[2] = -(p[2] + p[5] + p[7] + p[19])*u[2]            + 
              p[7]*u[3] + p[19]*u[8] + 2.0*p[2]*u[2]*u[5]
    # area AB
    du[3] = -(p[1] + p[2] + p[3] + p[8] + p[9] + p[19])*u[3] + 
              p[8]*u[1] + p[9]*u[2]+ p[19]*u[9]              +
              p[1]*(u[4]*u[3] + u[6]*u[1])                   + 
              p[2]*(u[5]*u[3] + u[6]*u[2])                   + 
              p[3]*(u[4]*u[2] + u[5]*u[1])
    ## extinction
    # area A
    du[4] = -(p[1] + p[4] + p[6] + p[19])*u[4]             + 
              p[4] + p[19]*u[10] + p[6]*u[6] + p[1]*u[4]^2
    # area B
    du[5] = -(p[2] + p[5] + p[7] + p[19])*u[5]             + 
              p[5] + p[19]*u[11] + p[7]*u[6] + p[2]*u[5]^2
    # area AB
    du[6] = -(p[1] + p[2] + p[3] + p[8] + p[9] + p[19])*u[6]   +
              p[8]*u[4] + p[9]*u[5] + p[19]*u[12]              + 
              p[1]*u[6]*u[4] + p[2]*u[6]*u[5] + p[3]*u[4]*u[5]
    ### Hidden State 2
    ## probabilities 
    # area A
    du[7] = -(p[10] + p[13] + p[15] + p[20])*u[7]            + 
              p[15]*u[9] + p[20]*u[1] + 2.0*p[10]*u[7]*u[10]
    # area B
    du[8] = -(p[11] + p[14] + p[16] + p[20])*u[8]            + 
              p[16]*u[9] + p[20]*u[2] + 2.0*p[11]*u[8]*u[11]
    # area AB
    du[9] = -(p[10] + p[11] + p[12] + p[17] + p[18] + p[20])*u[9] + 
              p[17]*u[7] + p[18]*u[8]+ p[20]*u[3]                 +
              p[10]*(u[10]*u[9] + u[12]*u[7])                     + 
              p[11]*(u[11]*u[9] + u[12]*u[8])                     + 
              p[12]*(u[10]*u[8] + u[11]*u[7])
    ## extinction
    # area A
    du[10] = -(p[10] + p[13] + p[15] + p[20])*u[10]            + 
               p[13] + p[20]*u[4] + p[15]*u[12] + p[10]*u[10]^2
    # area B
    du[11] = -(p[11] + p[14] + p[16] + p[20])*u[11]            + 
               p[14] + p[20]*u[5] + p[16]*u[12] + p[11]*u[11]^2
    # area AB
    du[12] = -(p[10] + p[11] + p[12] + p[17] + p[18] + p[20])*u[6]       +
               p[17]*u[10] + p[18]*u[11] + p[20]*u[6]                    + 
               p[10]*u[10]*u[12] + p[11]*u[11]*u[12] + p[12]*u[10]*u[11]
  end

  return nothing
end





#=
  Model with separate endemic extinction rates.

      1    2    3     4    5    6         7         8         9       
  p = sa1, sb1, sab1, ma1, mb1, qa1->ab1, qb1->ab1, qab1->a1, qab1->b1,

      10   11   12    13   14   15        16        17        18      
      sa2, sb2, sab2, ma2, mb2, qa2->ab2, qb2->ab2, qab2->a2, qab2->b2,

      19       20       21   22   23   24
      wi1->i2, wi2->i1, βa1, βb1, βa2, βb2,

      1    2    3     4    5    6     7    8    9     10   11   12
  u = da1, db1, dab1, ea1, eb1, eab1, da2, db2, dab2, ea2, eb2, eab2 
=#
"""
    egeohisse_2k_s(du::Array{Float64,1}, 
                   u::Array{Float64,1}, 
                   p::Array{Float64,1}, 
                   t::Float64)

Speciation EGeoHiSSE equation for 2 areas and 2 hidden states.
"""
function make_egeohisse_2k_s(af::Function)

  function f(du::Array{Float64,1}, 
             u::Array{Float64,1}, 
             p::Array{Float64,1}, 
             t::Float64)

    # UNIdimensional time function
    aft = af(t)

    sa1 = p[1]  * exp(aft*p[21])
    sb1 = p[2]  * exp(aft*p[22])
    sa2 = p[10] * exp(aft*p[23])
    sb2 = p[11] * exp(aft*p[24])

    @inbounds begin

      ### Hidden State 1
      ## probabilities 
      # area A
      du[1] = -(sa1 + p[4] + p[6] + p[19])*u[1]            + 
                p[6]*u[3] + p[19]*u[7] + 2.0*sa1*u[1]*u[4]
      # area B
      du[2] = -(sb1 + p[5] + p[7] + p[19])*u[2]            + 
                p[7]*u[3] + p[19]*u[8] + 2.0*sb1*u[2]*u[5]
      # area AB
      du[3] = -(sa1 + sb1 + p[3] + p[8] + p[9] + p[19])*u[3] + 
                p[8]*u[1] + p[9]*u[2]+ p[19]*u[9]            +
                sa1*(u[4]*u[3] + u[6]*u[1])                  + 
                sb1*(u[5]*u[3] + u[6]*u[2])                  + 
                p[3]*(u[4]*u[2] + u[5]*u[1])
      ## extinction
      # area A
      du[4] = -(sa1 + p[4] + p[6] + p[19])*u[4]             + 
                p[4] + p[19]*u[10] + p[6]*u[6] + sa1*u[4]^2
      # area B
      du[5] = -(sb1 + p[5] + p[7] + p[19])*u[5]             + 
                p[5] + p[19]*u[11] + p[7]*u[6] + sb1*u[5]^2
      # area AB
      du[6] = -(sa1 + sb1 + p[3] + p[8] + p[9] + p[19])*u[6]   +
                p[8]*u[4] + p[9]*u[5] + p[19]*u[12]              + 
                sa1*u[6]*u[4] + sb1*u[6]*u[5] + p[3]*u[4]*u[5]
      ### Hidden State 2
      ## probabilities 
      # area A
      du[7] = -(sa2 + p[13] + p[15] + p[20])*u[7]            + 
                p[15]*u[9] + p[20]*u[1] + 2.0*sa2*u[7]*u[10]
      # area B
      du[8] = -(sb2 + p[14] + p[16] + p[20])*u[8]            + 
                p[16]*u[9] + p[20]*u[2] + 2.0*sb2*u[8]*u[11]
      # area AB
      du[9] = -(sa2 + sb2 + p[12] + p[17] + p[18] + p[20])*u[9] + 
                p[17]*u[7] + p[18]*u[8]+ p[20]*u[3]                 +
                sa2*(u[10]*u[9] + u[12]*u[7])                     + 
                sb2*(u[11]*u[9] + u[12]*u[8])                     + 
                p[12]*(u[10]*u[8] + u[11]*u[7])
      ## extinction
      # area A
      du[10] = -(sa2 + p[13] + p[15] + p[20])*u[10]            + 
                 p[13] + p[20]*u[4] + p[15]*u[12] + sa2*u[10]^2
      # area B
      du[11] = -(sb2 + p[14] + p[16] + p[20])*u[11]            + 
                 p[14] + p[20]*u[5] + p[16]*u[12] + sb2*u[11]^2
      # area AB
      du[12] = -(sa2 + sb2 + p[12] + p[17] + p[18] + p[20])*u[6]       +
                 p[17]*u[10] + p[18]*u[11] + p[20]*u[6]                    + 
                 sa2*u[10]*u[12] + sb2*u[11]*u[12] + p[12]*u[10]*u[11]
    end

    return nothing
  end

  return f
end





#=
  Model with separate endemic extinction rates.

      1    2    3     4    5    6         7         8         9       
  p = sa1, sb1, sab1, ma1, mb1, qa1->ab1, qb1->ab1, qab1->a1, qab1->b1,

      10   11   12    13   14   15        16        17        18      
      sa2, sb2, sab2, ma2, mb2, qa2->ab2, qb2->ab2, qab2->a2, qab2->b2,

      19       20       21   22   23   24
      wi1->i2, wi2->i1, βa1, βb1, βa2, βb2,

      1    2    3     4    5    6     7    8    9     10   11   12
  u = da1, db1, dab1, ea1, eb1, eab1, da2, db2, dab2, ea2, eb2, eab2 
=#
"""
    egeohisse_2k_e(du::Array{Float64,1}, 
                   u::Array{Float64,1}, 
                   p::Array{Float64,1}, 
                   t::Float64)

Speciation EGeoHiSSE equation for 2 areas and 2 hidden states.
"""
function make_egeohisse_2k_e(af::Function)

  function f(du::Array{Float64,1}, 
             u::Array{Float64,1}, 
             p::Array{Float64,1}, 
             t::Float64)

    # UNIdimensional time function
    aft = af(t)

    ma1 = p[4]  * exp(aft*p[21])
    mb1 = p[5]  * exp(aft*p[22])
    ma2 = p[13] * exp(aft*p[23])
    mb2 = p[14] * exp(aft*p[24])

    @inbounds begin

      ### Hidden State 1
      ## probabilities 
      # area A
      du[1] = -(p[1] + ma1 + p[6] + p[19])*u[1]            + 
                p[6]*u[3] + p[19]*u[7] + 2.0*p[1]*u[1]*u[4]
      # area B
      du[2] = -(p[2] + mb1 + p[7] + p[19])*u[2]            + 
                p[7]*u[3] + p[19]*u[8] + 2.0*p[2]*u[2]*u[5]
      # area AB
      du[3] = -(p[1] + p[2] + p[3] + p[8] + p[9] + p[19])*u[3] + 
                p[8]*u[1] + p[9]*u[2]+ p[19]*u[9]              +
                p[1]*(u[4]*u[3] + u[6]*u[1])                   + 
                p[2]*(u[5]*u[3] + u[6]*u[2])                   + 
                p[3]*(u[4]*u[2] + u[5]*u[1])
      ## extinction
      # area A
      du[4] = -(p[1] + ma1 + p[6] + p[19])*u[4]             + 
                ma1 + p[19]*u[10] + p[6]*u[6] + p[1]*u[4]^2
      # area B
      du[5] = -(p[2] + mb1 + p[7] + p[19])*u[5]             + 
                mb1 + p[19]*u[11] + p[7]*u[6] + p[2]*u[5]^2
      # area AB
      du[6] = -(p[1] + p[2] + p[3] + p[8] + p[9] + p[19])*u[6]   +
                p[8]*u[4] + p[9]*u[5] + p[19]*u[12]              + 
                p[1]*u[6]*u[4] + p[2]*u[6]*u[5] + p[3]*u[4]*u[5]
      ### Hidden State 2
      ## probabilities 
      # area A
      du[7] = -(p[10] + ma2 + p[15] + p[20])*u[7]            + 
                p[15]*u[9] + p[20]*u[1] + 2.0*p[10]*u[7]*u[10]
      # area B
      du[8] = -(p[11] + mb2 + p[16] + p[20])*u[8]            + 
                p[16]*u[9] + p[20]*u[2] + 2.0*p[11]*u[8]*u[11]
      # area AB
      du[9] = -(p[10] + p[11] + p[12] + p[17] + p[18] + p[20])*u[9] + 
                p[17]*u[7] + p[18]*u[8]+ p[20]*u[3]                 +
                p[10]*(u[10]*u[9] + u[12]*u[7])                     + 
                p[11]*(u[11]*u[9] + u[12]*u[8])                     + 
                p[12]*(u[10]*u[8] + u[11]*u[7])
      ## extinction
      # area A
      du[10] = -(p[10] + ma2 + p[15] + p[20])*u[10]            + 
                 ma2 + p[20]*u[4] + p[15]*u[12] + p[10]*u[10]^2
      # area B
      du[11] = -(p[11] + mb2 + p[16] + p[20])*u[11]            + 
                 mb2 + p[20]*u[5] + p[16]*u[12] + p[11]*u[11]^2
      # area AB
      du[12] = -(p[10] + p[11] + p[12] + p[17] + p[18] + p[20])*u[6]       +
                 p[17]*u[10] + p[18]*u[11] + p[20]*u[6]                    + 
                 p[10]*u[10]*u[12] + p[11]*u[11]*u[12] + p[12]*u[10]*u[11]
    end

    return nothing
  end

  return f
end





