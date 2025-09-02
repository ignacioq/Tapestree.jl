using Test
using Tapestree
using Random
using DelimitedFiles

# read tree
tree = read_newick(joinpath(dirname(pathof(Tapestree)), "..", "data", "tree_5.tre"))
@test isa(tree, sT_label)
@test ntips(tree) === 5

#=
Constant Pure-birth
=# 
# simulate 
tr = sim_cpb(1.0, 0.5)
@test isa(tr, sTpb)

# perform inference
r, tv = insane_cpb(tree,
                   nburn  = 2,
                   niter  = 5,
                   nthin  = 5,
                   nflush = 5,
                   ofile  = homedir()*"/test")

@test isa(r, Matrix{Float64})
@test isa(tv, Vector{sTpb})


#=
Constant Birth-Death
=# 
# simulate 
tr = sim_cbd(1.0, 1.0, 0.3)
@test isa(tr, sTbd)

# perform inference
r, tv = insane_cbd(tree,
                   nburn  = 2,
                   niter  = 5,
                   nthin  = 5,
                   nflush = 5,
                   ofile  = homedir()*"/test")

@test isa(r, Matrix{Float64})
@test isa(tv, Vector{sTbd})


#=
Constant Fossilized Birth-Death
=# 
# read tree
ftree = read_newick(joinpath(dirname(pathof(Tapestree)), "..", "data", "tree_6.tre"), true)
@test isa(ftree, sTf_label)
@test ntips(ftree) === 8


# simulate 
tr = sim_cfbd(1.0, 1.0, 0.3, 0.3)
@test isa(tr, sTfbd)

# perform inference
r, tv = insane_cfbd(ftree,
                    nburn  = 2,
                    niter  = 5,
                    nthin  = 5,
                    nflush = 5,
                    ofile  = homedir()*"/test")

@test isa(r, Matrix{Float64})
@test isa(tv, Vector{sTfbd})


#=
Constant Occurrence Birth-Death
=# 
# read tree
ftree = read_newick(joinpath(dirname(pathof(Tapestree)), "..", "data", "tree_6.tre"), true)
ωtimes = readdlm(joinpath(dirname(pathof(Tapestree)), "..", "data", "fossil_occurrences.csv"), ';')[:]
@test isa(ftree, sTf_label)
@test isa(ωtimes, Vector{Float64})
@test ntips(ftree) === 8
@test length(ωtimes) === 50

# simulate 
tr, occ = sim_cobd(1.0, 1.0, 0.3, 0.3, 5.0)
@test isa(tr, sTfbd)
@test isa(occ, Vector{Float64})

# perform inference
r, tv = insane_cobd(ftree,
                    ωtimes,
                    nburn   = 2,
                    niter   = 5,
                    nthin   = 5,
                    nflushθ = 5,
                    nflushΞ = 5,
                    ofile   = homedir()*"/test")

@test isa(r, Matrix{Float64})
@test isa(tv, Vector{sTfbd})


#=
Birth-Death Diffusion: Pure-birth
=# 

# simulate based on number of species
Random.seed!(7)
tr = sim_gbmpb(20, λ0 = .5, α = 0.0, σλ = 0.1)
@test isa(tr, iTpb)
@test ntips(tr) === 20

# simulate based on time
Random.seed!()
tr = sim_gbmpb(10.0, λ0 = .5, α = 0.0, σλ = 0.1)
@test isa(tr, iTpb)

# perform inference
r, tv = insane_gbmpb(tree,
                     nburn  = 2,
                     niter  = 5,
                     nthin  = 5,
                     nflush = 5,
                     ofile  = homedir()*"/test")

@test isa(r, Matrix{Float64})
@test isa(tv, Vector{iTpb})



#=
Birth-Death Diffusion: constant death
=# 

# simulate based on time
tr = sim_gbmce(10.0, λ0 = .5, α = 0.0, μ = 0.1, σλ = 0.1)
@test isa(tr, iTce)

# perform inference
r, tv = insane_gbmce(tree,
                     nburn  = 2,
                     niter  = 5,
                     nthin  = 5,
                     nflush = 5,
                     ofile  = homedir()*"/test")

@test isa(r, Matrix{Float64})
@test isa(tv, Vector{iTce})



#=
Birth-Death Diffusion: constant turnover
=# 

# simulate based on time
tr = sim_gbmct(10.0, λ0 = .5, α = 0.0, ϵ = 0.1, σλ = 0.1)
@test isa(tr, iTct)

# perform inference
r, tv = insane_gbmct(tree,
                     nburn  = 2,
                     niter  = 5,
                     nthin  = 5,
                     nflush = 5,
                     ofile  = homedir()*"/test")

@test isa(r, Matrix{Float64})
@test isa(tv, Vector{iTct})


#=
Birth-Death Diffusion
=# 

# simulate based on time
tr = sim_gbmbd(10.0, λ0 = .5, μ0 = .2, α = 0.0, σλ = 0.1, σμ = 0.1)
@test isa(tr, iTbd)

# perform inference
r, tv = insane_gbmbd(tree,
                     nburn  = 2,
                     niter  = 5,
                     nthin  = 5,
                     nflush = 5,
                     ofile  = homedir()*"/test")

@test isa(r, Matrix{Float64})
@test isa(tv, Vector{iTbd})



#=
Fossilized Birth-Death Diffusion
=# 

# simulate based on time
tr = sim_gbmfbd(10.0, λ0 = .5, μ0 = .2, αλ = 0.0, αμ = 0.0, σλ = 0.1, σμ = 0.1)
@test isa(tr, iTfbd)

# perform inference
r, tv = insane_gbmfbd(ftree,
                      nburn  = 2,
                      niter  = 5,
                      nthin  = 5,
                      nflush = 5,
                      ofile  = homedir()*"/test")

@test isa(r, Matrix{Float64})
@test isa(tv, Vector{iTfbd})



#=
Fossilized Birth-Death Diffusion
=# 

# simulate based on time
tr, occ = sim_gbmobd(10.0, λ0 = .5, μ0 = .2, αλ = 0.0, αμ = 0.0, σλ = 0.1, σμ = 0.1)
@test isa(tr, iTfbd)
@test isa(occ, Vector{Float64})

# perform inference
r, tv = insane_gbmobd(ftree,
                      ωtimes,
                      nburn   = 2,
                      niter   = 5,
                      nthin   = 5,
                      nflushθ = 5,
                      nflushΞ = 5,
                      ofile  = homedir()*"/test")

@test isa(r, Matrix{Float64})
@test isa(tv, Vector{iTfbd})



#=
Diffused Brownian motion (DBM)
=# 

# simulate based on time
tr, xav = sim_dbm(ftree, 0.0, 0.0, 0.1, 0.0, 0.1, 1e-3)
@test isa(tr, sTxs)
@test isa(xav, Dict{String, Float64})
@test length(xav) === 10

# perform inference
r, tv = insane_dbm(ftree, xav,
                   nburn  = 2,
                   niter  = 5,
                   nthin  = 5,
                   nflush = 5,
                   ofile  = homedir()*"/test")

@test isa(r, Matrix{Float64})
@test isa(tv, Vector{sTxs})

rm(homedir()*"/test.txt") 
rm(homedir()*"/test.log") 


# ltt
ltti = ltt(tree)

@test isa(ltti, Ltt)
@test isa(ltti.n, Vector{Int64})

ltti = ltt(ftree)

@test isa(ltti, Ltt)
@test isa(ltti.t, Vector{Float64})





