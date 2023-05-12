var documenterSearchIndex = {"docs":
[{"location":"man/esse/#ESSE","page":"ESSE","title":"ESSE","text":"","category":"section"},{"location":"man/esse/#Example","page":"ESSE","title":"Example","text":"","category":"section"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"In this example we run some random data for 50 species across 2 areas, where each area has a specific covariate that affects speciation and we also allow for 2 hidden states.","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"Open Julia and load the Tapestree package: ","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"using Tapestree","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"Specify the path to the phylogenetic tree (in a format that ape::read.tree() can read):","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"tree_file    = joinpath(dirname(pathof(Tapestree)), \"..\", \"data\", \"tree_50.tre\")","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"Specify state data. Data should be a .txt file where each row is a species,  first the species name that matches the tree tip labels and each subsequent column specifying presence or absence in a given area with 0 or 1,  respecitvely. Open st2_data.txt in the data folder to see an example for 2 areas.","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"states_file = joinpath(dirname(pathof(Tapestree)), \"..\", \"data\", \"st2_data.txt\")","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"Specify covariate data y = f(x). Data should be a .txt file where the  first column is time x in backward fashion (the present is 0 and the  past is  0), and the subsequent columns are the respective time covariates for each area f(x). If there is only one covariate, the same is used across all areas, if not, the number of covariates should match the specific model.  More than one covariate per area is allowed, and in the case of covariates  affecting colonization rates, they should match the number of possible  colonization parameters between all areas. Open env_data_2.txt in the data folder to see an example for 2 covariates for 2 areas.","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"states_file = joinpath(dirname(pathof(Tapestree)), \"..\", \"data\", \"st2_data.txt\")","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"Specify output MCMC file (homedir() is an alias to your home folder)","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"out_file  = *(homedir(),\"...\")","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"Specify the (optional) output file (homedir() is an alias to your home  folder) for the node marginal probabilities. ","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"out_states  = *(homedir(),\"...\")","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"Run the esse() (ESSE: Environmental and State dependent Speciation and  Extinction) model, with covariates affecting speciation rates:","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"esse(tree_file, out_file, 2, envdata_file = envdata_file, \n  states_file = states_file, out_states = out_states, cov_mod = (\"s\",))","category":"page"},{"location":"man/esse/#Parallel-MC3","page":"ESSE","title":"Parallel MC3","text":"","category":"section"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"It is encouraged to use Metropolis coupled MCMC (MC3) for more robust  convergence (the posterior surfaced is highly peaked).","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"Load the Distributed package, set the number of processors for Julia, and make Tapestree available to all (see the Distributed  package for more information). Below we add 3 processors.","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"using Distributed\naddprocs(3)\n@everywhere using Tapestree","category":"page"},{"location":"man/esse/#Set-parameter-constraints","page":"ESSE","title":"Set parameter constraints","text":"","category":"section"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"To constrain parameters to be equal to one another or to fix them to be 0,  it is necessary to create equalities between parameters and pass them as an  argument. In general parameters are specified as follows:    \"<parameter name>_<area>_<hidden state>\"","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"Speciation is \"lambda\" (e.g. speciation rate for area A and hidden state 0: lambda_A_0)\nLocal extinction is \"loss\" (e.g. local extinction rate for area B and hidden state 1: loss_B_1)\ncolonization is \"gain\" (e.g. gain rate from A -> B and hidden state 1: gain_AB_1)\nthe effect of the covariate is \"beta_<effect parameter name>\" (e.g. effect of first covariate on speciation in area A and hidden state 0: beta_lambda_1_A_0 - here the 1 after the lambda is because is the first covariate)\ntransition between hidden states is \"q\" (e.g. transition rate from hidden state 0 -> 1 and hidden state 1: q_01). ","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"You have to make sure that the given constraints apply to the specification of  model. In the following example, we constraint local and global extinction rates to be the same across 2 areas and 2 hidden states , and constrain the hidden states transition rates.","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"cpar = (\"q_01 = q_10\",\n        \"loss_A_0 = mu_A_0\",\n        \"loss_B_0 = mu_B_0\",\n        \"loss_A_1 = mu_A_1\",\n        \"loss_B_1 = mu_B_1\")","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"We now run esse, using MC3 with Metropolis-Hastings MCMC (here using 3 parallel chains).","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"esse(tree_file, out_file, 2,\n     envdata_file = envdata_file,\n     states_file  = states_file, \n     out_states   = out_states,\n     constraints  = cpar,\n     cov_mod      = (\"s\",),\n     ncch         = 3,\n     parallel     = true,\n     niter        = 5_000,\n     nthin        = 100,\n     dt           = 0.8,\n     nburn        = 1_000, \n     mc           = \"mh\",\n     node_ps      = (true, 100))","category":"page"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"For optional (keyword) arguments, see below","category":"page"},{"location":"man/esse/#Function-Documentation","page":"ESSE","title":"Function Documentation","text":"","category":"section"},{"location":"man/esse/","page":"ESSE","title":"ESSE","text":"esse\nsimulate_sse","category":"page"},{"location":"man/esse/#Tapestree.ESSE.esse","page":"ESSE","title":"Tapestree.ESSE.esse","text":"esse(tree_file   ::String,\n     out_file    ::String,\n     h           ::Int64;\n     states_file ::String            = \"NaN\",\n     envdata_file::String            = \"NaN\",\n     cov_mod     ::NTuple{M,String}  = (\"\",),\n     node_ps     ::Tuple{Bool,Int64} = (true, 10),\n     out_states  ::String            = \"\",\n     constraints ::NTuple{N,String}  = (\" \",),\n     mvpars      ::NTuple{O,String}  = (\" \",),\n     niter       ::Int64             = 10_000,\n     nthin       ::Int64             = 10,\n     nburn       ::Int64             = 200,\n     tune_int    ::Int64             = 100,\n     nswap       ::Int64             = 10,\n     ncch        ::Int64             = 1,\n     parallel    ::Bool              = ncch > 1,\n     dt           ::Float64          = 0.2,\n     ntakew      ::Int64             = 100,\n     winit       ::Float64           = 2.0,\n     scale_y     ::NTuple{2,Bool}    = (true, false),\n     algorithm   ::String            = \"pruning\",\n     mc          ::String            = \"slice\",\n     λpriors     ::Float64           = .1,\n     μpriors     ::Float64           = .1,\n     gpriors     ::Float64           = .1,\n     lpriors     ::Float64           = .1,\n     qpriors     ::Float64           = .1,\n     βpriors     ::NTuple{2,Float64} = (0.0, 10.0),\n     hpriors     ::Float64           = .1,\n     optimal_w   ::Float64           = 0.8,\n     tni         ::Float64           = 1.0,\n     obj_ar      ::Float64           = 0.6,\n     screen_print::Int64             = 5,\n     Eδt         ::Float64           = 1e-3,\n     ti          ::Float64           = 0.0,\n     ρ           ::Array{Float64,1}  = [1.0]) where {M,N,O}\n\nRun geographic esse. See tutorial for how these files should be specified.\n\n...\n\nArguments\n\ntree_file::String: full path to tree file.\nout_file::String: full path to write MCMC output.\nh::Int64: number of hidden states.\nstates_file ::String = \"NaN\": full path to states file. If \"NaN\", no \n\nobserved states are used (only hidden states).\n\nenvdata_file::String = \"NaN\": full path to covariates file. If \"NaN\", no \n\ncovariates are used (i.e., constant rates).\n\ncov_mod::NTuple{M,String} = (\"\",): specifies which rates are affected by \n\ncovariates: s for speciation, e for extinction, and g for colonization. More than 1 is possible.\n\nnode_ps::Tuple{Bool,Int64} = (true, 10): first index specifies if posterior\n\nmarginal probabilities for nodes should be computed, second index the number of iterations to be computed. \n\nout_states::String = \"\": full path to write node probabilities output.\nconstraints::NTuple{N,String} = (\" \",): constraints for the model \n\nparameters.\n\nmvpars::NTuple{O,String} = (\" \",): which parameters should be multivariate \n\nwhen using slice sampling for better convergence.\n\nniter::Int64 = 10_000: number of iterations.\nnthin::Int64 = 10: frequency at which to record MCMC state.\nnburn::Int64 = 200: number of iterations to discard as burn-in.\ntune_int::Int64 = 100: number of iterations during nburn to tune proposal\n\nwindow for MH.\n\nnswap::Int64 = 10: every iteration to try to swap chain likelihoods in MC3.\nncch::Int64 = 1: number of chains.\nparallel::Bool = false: if parallel run.\ndt::Float64 = 0.2: temperature for MC3.\nntakew::Int64 = 100: number of iterations from nburn to tune the window\n\nfor slice sampling.\n\nwinit::Float64 = 2.0: initial window for slice sampling.\nscale_y::NTuple{2,Bool} = (true, false): first index if scale covariates y\n\nto [0,1], second, if scale covariates y all together between [0,1].\n\nalgorithm::String = \"pruning\": likelihood algorithm between pruning \n\n(recommended) or flow.\n\nmc::String = \"slice\": which sampling slice (slice-sampling) or \n\nmh (metropolis-hasting).\n\nλpriors::Float64 = 0.1: rate of Exponential prior for speciation.\nμpriors::Float64 = 0.1: rate of Exponential prior for global extinction.\ngpriors::Float64 = 0.1: rate of Exponential prior for colonization.\nlpriors::Float64 = 0.1: rate of Exponential prior for local extinction.\nqpriors::Float64 = 0.1: rate of Exponential prior for hidden state \n\ntransitions.\n\nβpriors::NTuple{2,Float64} = (0.0, 10.0): mean and variance of Normal \n\nprior for effect of covariates.\n\nhpriors::Float64 = 0.1: rate of Exponential prior for differences between\n\nhidden states.\n\noptimal_w::Float64 = 0.8: optimal window.\ntni::Float64 = 1.0: initial tuning for rates.\nobj_ar::Float64 = 0.23: objective acceptance rate.\nscreen_print::Int64 = 5: seconds to wait to update screen log.\nEδt::Float64 = 1e-3: for flow algorithm.\nti::Float64 = 0.0: for flow algorithm.\nρ::Array{Float64,1} = [1.0]: sampling fraction for each state (each area and\n\nwidespread). ...\n\n...\n\nReturned values\n\nArray of the mcmc parameters.\n\n...\n\n\n\n\n\nesse(tv          ::Dict{Int64,Array{Float64,1}},\n     ed          ::Array{Int64,2}, \n     el          ::Array{Float64,1}, \n     x           ::Array{Float64,1},\n     y           ::Array{Float64,L}, \n     cov_mod     ::NTuple{M,String},\n     out_file    ::String,\n     h           ::Int64;\n     constraints ::NTuple{N,String}  = (\" \",),\n     mvpars      ::NTuple{O,String}  = (\"lambda = beta\",),\n     niter       ::Int64             = 10_000,\n     nthin       ::Int64             = 10,\n     nburn       ::Int64             = 200,\n     ncch        ::Int64             = 1,\n     ntakew      ::Int64             = 100,\n     winit       ::Float64             = 2.0,\n     scale_y     ::NTuple{2,Bool}    = (true, false),\n     algorithm   ::String            = \"pruning\",\n     λpriors     ::Float64           = .1,\n     μpriors     ::Float64           = .1,\n     gpriors     ::Float64           = .1,\n     lpriors     ::Float64           = .1,\n     qpriors     ::Float64           = .1,\n     βpriors     ::NTuple{2,Float64} = (0.0, 10.0),\n     hpriors     ::Float64           = .1,\n     optimal_w   ::Float64           = 0.8,\n     screen_print::Int64             = 5,\n     Eδt         ::Float64           = 1e-3,\n     ti          ::Float64           = 0.0,\n     ρ           ::Array{Float64,1}  = [1.0]) where {L,M,N,O}\n\nWrapper for running a esse model from simulations.\n\n\n\n\n\n","category":"function"},{"location":"man/esse/#Tapestree.ESSE.simulate_sse","page":"ESSE","title":"Tapestree.ESSE.simulate_sse","text":"simulate_sse(λ          ::Array{Float64,1},\n             μ          ::Array{Float64,1},\n             l          ::Array{Float64,1},\n             g          ::Array{Float64,1},\n             q          ::Array{Float64,1},\n             t          ::Float64;\n             δt         ::Float64 = 1e-4,\n             ast        ::Int64   = 0,\n             nspp_max   ::Int64   = 100_000,\n             retry_ext  ::Bool    = true,\n             rejectel0  ::Bool    = true,\n             verbose    ::Bool    = true,\n             rm_ext     ::Bool    = true,\n             states_only::Bool    = false, \n             start      ::Symbol  = :crown)\n\nSimulate tree according to the geographic esse model. The number of areas  and hidden states is inferred from the parameter vectors, but they must be  consistent between them and with the covariates. See tutorial for an example.\n\n...\n\nArguments\n\nλ::Array{Float64,1}: rates for within-area and between-area speciation.\nμ::Array{Float64,1}: per-area extinction rates when it leads to global \n\nextinction.\n\nl::Array{Float64,1}: per-area extinction rates when it leads to local \n\nextinction.\n\ng::Array{Float64,1}: colonization rates between areas. \nq::Array{Float64,1}: transition rates between hidden states.\nt::Float64: simulation time.\nδt::Float64 = 1e-4: time step to perform simulations. Smaller more precise \n\nbut more computationally expensive.\n\nast::Int64 = 0: initial state. 0 specifies random sampling based on the \n\ninput parameters.\n\nnspp_max::Int64 = 100_000: maximum number of species allowed to stop \n\nsimulation.\n\nretry_ext::Bool = true: automatically restart simulation if simulation goes\n\nextinct.\n\nrejectel0::Bool = true: reject simulations where there are edges of 0 \n\nlength.\n\nverbose::Bool = true: print messages.\nrm_ext::Bool = true: remove extinct taxa from output.\nstates_only::Bool = false: if only return tip states (faster).\nstart::Symbol  = :crown: if crown, starts after a speciation event with \n\ntwo lineages, if stem, starts with one lineage. ...\n\n...\n\nReturned values\n\nDictionary with tip number and corresponding state.\nArray with parent -> daughter edges.\nArray with edge lengths.\nNumber of maximum species.\n\n...\n\n\n\n\n\nsimulate_sse(λ          ::Array{Float64,1},\n             μ          ::Array{Float64,1},\n             l          ::Array{Float64,1},\n             g          ::Array{Float64,1},\n             q          ::Array{Float64,1},\n             β          ::Array{Float64,1},\n             t          ::Float64,\n             x          ::Array{Float64,1},\n             y          ::Array{Float64,N},\n             cov_mod    ::Tuple{Vararg{String}};\n             δt         ::Float64 = 1e-4,\n             ast        ::Int64   = 0,\n             nspp_max   ::Int64   = 100_000,\n             retry_ext  ::Bool    = true,\n             rejectel0  ::Bool    = true,\n             verbose    ::Bool    = true,\n             rm_ext     ::Bool    = true,\n             states_only::Bool    = false,\n             start      ::Symbol  = :crown) where {N}\n\nSimulate tree according to the geographic sse model. The number of areas  and hidden states is inferred from the parameter vectors, but they must be  consistent between them and with the covariates. See tutorial for an example.\n\n...\n\nArguments\n\nλ::Array{Float64,1}: rates for within-area and between-area speciation.\nμ::Array{Float64,1}: per-area extinction rates when it leads to global \n\nextinction.\n\nl::Array{Float64,1}: per-area extinction rates when it leads to local \n\nextinction.\n\ng::Array{Float64,1}: colonization rates between areas. \nq::Array{Float64,1}: transition rates between hidden states.\nβ::Array{Float64,1}: per-area effect of covariates.\nt::Float64: simulation time.\nx::Array{Float64,1}: times where the covariate y is sampled.\ny::Array{Float64,N}: value of covariates, i.e., f(x). Can be multivariate.\ncov_mod::Tuple{Vararg{String}}: specifies which rates are affected by \n\ncovariates: s for speciation, e for extinction, and g for colonization. More than 1 is possible. \n\nδt::Float64 = 1e-4: time step to perform simulations. Smaller more precise \n\nbut more computationally expensive.\n\nast::Int64 = 0: initial state. 0 specifies random sampling based on the \n\ninput parameters.\n\nnspp_max::Int64 = 100_000: maximum number of species allowed to stop \n\nsimulation.\n\nretry_ext::Bool = true: automatically restart simulation if simulation goes\n\nextinct.\n\nrejectel0::Bool = true: reject simulations where there are edges of 0 \n\nlength.\n\nverbose::Bool = true: print messages.\nrm_ext::Bool = true: remove extinct taxa from output.\nstates_only::Bool = false: if only return tip states (faster).\nstart::Symbol  = :crown: if crown, starts after a speciation event with \n\ntwo lineages, if stem, starts with one lineage. ...\n\n...\n\nReturned values\n\nDictionary with tip number and corresponding state.\nArray with parent -> daughter edges.\nArray with edge lengths.\nNumber of maximum species.\n\n...\n\n\n\n\n\n","category":"function"},{"location":"man/tribe/#TRIBE","page":"TRIBE","title":"TRIBE","text":"","category":"section"},{"location":"man/tribe/#Example","page":"TRIBE","title":"Example","text":"","category":"section"},{"location":"man/tribe/","page":"TRIBE","title":"TRIBE","text":"Open Julia and load the Tapestree package: ","category":"page"},{"location":"man/tribe/","page":"TRIBE","title":"TRIBE","text":"using Tapestree","category":"page"},{"location":"man/tribe/","page":"TRIBE","title":"TRIBE","text":"Specify the path to the phylogenetic tree (in a format that ape::read.tree() can read):","category":"page"},{"location":"man/tribe/","page":"TRIBE","title":"TRIBE","text":"finches_tree_file = joinpath(dirname(pathof(Tapestree)), \"..\", \"data\", \"finches_rescaled.tre\")","category":"page"},{"location":"man/tribe/","page":"TRIBE","title":"TRIBE","text":"Specify data. Data should be a .txt file where each row is a species, first  the species name that matches the tree tip labels, second the phenotypic data  and then the species presence in each area (0 if absent and 1 if present).  Open finches_pca1.txt in the data folder to see an example.","category":"page"},{"location":"man/tribe/","page":"TRIBE","title":"TRIBE","text":"finches_data_file = joinpath(dirname(pathof(Tapestree)), \"..\", \"data\", \"finches_pca1.txt\")","category":"page"},{"location":"man/tribe/","page":"TRIBE","title":"TRIBE","text":"Specify output file (homedir() is an alias to your home folder)","category":"page"},{"location":"man/tribe/","page":"TRIBE","title":"TRIBE","text":"out_file  = *(homedir(),\"...\")","category":"page"},{"location":"man/tribe/","page":"TRIBE","title":"TRIBE","text":"Run the tribe() (TRIBE: Trait and Range Interspecific Biogeographic Evolution) model:","category":"page"},{"location":"man/tribe/","page":"TRIBE","title":"TRIBE","text":"tribe(finches_tree_file, finches_data_file, out_file)","category":"page"},{"location":"man/tribe/","page":"TRIBE","title":"TRIBE","text":"For optional (keyword) arguments, see below","category":"page"},{"location":"man/tribe/#Function-Documentation","page":"TRIBE","title":"Function Documentation","text":"","category":"section"},{"location":"man/tribe/","page":"TRIBE","title":"TRIBE","text":"tribe\nsimulate_tribe","category":"page"},{"location":"man/tribe/#Tapestree.TRIBE.tribe","page":"TRIBE","title":"Tapestree.TRIBE.tribe","text":"tribe(tree_file   ::String,\n      data_file   ::String,\n      out_file    ::String;\n      min_dt      ::Float64           = 0.01,\n      niter       ::Int64             = 50_000,\n      nburn       ::Int64             = 50_000,\n      nthin       ::Int64             = 500,\n      saveXY      ::Tuple{Bool,Int64} = (false, 1_000),\n      saveDM      ::Tuple{Bool,Int64} = (false, 1_000),\n      ωxprior     ::NTuple{2,Float64} = (0.,10.),\n      ω1prior     ::NTuple{2,Float64} = (0.,10.),\n      ω0prior     ::NTuple{2,Float64} = (0.,10.),\n      σ²prior     ::Float64           = 1e-1,\n      λprior      ::Float64           = 1e-1,\n      weight      ::NTuple{5,Float64} = (0.15,0.05,0.02,0.02,5e-3),\n      λ1i         ::Float64           = 1.0,\n      λ0i         ::Float64           = 0.5,\n      ωxi         ::Float64           = 0.0,\n      ω1i         ::Float64           = 0.0,\n      ω0i         ::Float64           = 0.0,\n      fix_ωx      ::Bool              = false,\n      fix_ω1      ::Bool              = false,\n      fix_ω0      ::Bool              = false,\n      delim       ::Char              = '\t',\n      eol         ::Char              = '\r',\n      screen_print::Int64             = 5)\n\nRun tribe model. \n\n...\n\nArguments\n\ntree_file::String: full path to tree file.\ndata_file::String: full path to data file.\nout_file::String: full path to write MCMC output.\nmin_dt::Float64 = 0.01: the percentage of tree height allowed for \n\ndiscretization (lower values are more precise but take longer).\n\nniter::Int64 = 50_000: the number of iterations.\nnburn::Int64 = 50_000: the number of iterations in the adaptive burn-in phase.\nnthin::Int64 = 500: the iteration sampling frequency.\nsaveXY::Tuple{Bool,Int64} = (false, 1_000): first index to \n\nsave (or not) data augmented histories, second index for sampling frequency.\n\nsaveDM::Tuple{Bool,Int64} = (false, 1_000): a tuple of length 2: first is a boolean to save (or not) data augmented deterministic effects, second an integer for sampling frequency.\nωxprior::NTuple{2,Float64} = (0.,10.): a tuple of length 2 for the normal prior of ωx, first the mean, second the variance.\nω1prior::NTuple{2,Float64} = (0.,10.): a tuple of length 2 for the normal prior of ω1, first the mean, second the variance.\nω0prior::NTuple{2,Float64} = (0.,10.): a tuple of length 2 for the normal prior of ω0, first the mean, second the variance.\nσ²prior::Float64 = 1e-1: a float for the mean of the exponential prior for σ².\nλprior::Float64 = 1e-1: a float for the mean of the exponential prior for both λs.\nweight::NTuple{5,Float64} = (0.15,0.05,0.02,0.02,5e-3): a tuple of length 5 specifying the probabilities to update σ², ωx, ω1 & ω0, and λ1 & λ0 respectively.\nλ1i::Float64 = 1.0: a float for the starting value for λ1.\nλ0i::Float64 = 0.5: a float for the starting value for λ0.\nωxi::Float64 = 0.0: a float for the starting value for ωx.\nω1i::Float64 = 0.0: a float for the starting value for ω1.\nω0i::Float64 = 0.0: a float for the starting value for ω0.\nfix_ωx::Bool = false: a boolean to make inference without ωx.\nfix_ω1::Bool = false: a boolean to make inference without ω1.\nfix_ω0::Bool = false: a boolean to make inference without ω0.\ndelim::Char= '\t': for ddlm.\neol::Char= ' ': for ddlm.\nscreen_print::Int64 = 5: seconds to wait to update screen log.\n\n...\n\n...\n\nReturned values\n\nArray of the mcmc parameters.\n\n...\n\n\n\n\n\ntribe(tip_values::Dict{Int64,Float64}, \n      tip_areas ::Dict{Int64,Array{Int64,1}},\n      tree      ::rtree, \n      bts       ::Array{Float64,1},\n      out_file  ::String;\n      min_dt    ::Float64           = 0.01,\n      niter     ::Int64             = 500_000,\n      nburn     ::Int64             = 500_000,\n      nthin     ::Int64             = 1_000,\n      saveXY    ::Tuple{Bool,Int64} = (false, 1_000),\n      saveDM    ::Tuple{Bool,Int64} = (false, 1_000),\n      ωxprior   ::NTuple{2,Float64} = (0.,10.),\n      ω1prior   ::NTuple{2,Float64} = (0.,10.),\n      ω0prior   ::NTuple{2,Float64} = (0.,10.),\n      σ²prior   ::Float64           = 1e-1,\n      λprior    ::Float64           = 1e-1,\n      weight    ::NTuple{5,Float64} = (0.15,0.05,0.02,0.02,5e-3),\n      λ1i       ::Float64           = 1.0,\n      λ0i       ::Float64           = 0.4,\n      ωxi       ::Float64           = 0.,\n      ω1i       ::Float64           = 0.,\n      ω0i       ::Float64           = 0.,\n      fix_ωx    ::Bool              = false,\n      fix_ω1    ::Bool              = false,\n      fix_ω0    ::Bool              = false,\n      delim     ::Char              = '\t',\n      eol       ::Char              = '\r')\n\nRun tribe for simulations. Wrapper for all functions.\n\n\n\n\n\ntribe(out_file::String;\n      niter   ::Int64             = 500_000,\n      nburn   ::Int64             = 500_000,\n      nthin   ::Int64             = 1_000,\n      ωxprior ::NTuple{2,Float64} = (0.,10.),\n      ω1prior ::NTuple{2,Float64} = (0.,10.),\n      ω0prior ::NTuple{2,Float64} = (0.,10.),\n      σ²prior ::Float64           = 1e-1,\n      λprior  ::Float64           = 1e-1,\n      weight  ::NTuple{4,Float64} = (0.15,0.05,0.02,0.02),\n      σ²i     ::Float64           = 1.,\n      ωxi     ::Float64           = 0.,\n      ω1i     ::Float64           = 0.01,\n      ω0i     ::Float64           = 0.01,\n      λ1i     ::Float64           = 1.0,\n      λ0i     ::Float64           = 0.2,\n      fix_ωx  ::Bool              = false,\n      fix_ω1  ::Bool              = false,\n      fix_ω0  ::Bool              = false)\n\nRun tribe under the prior. Wrapper for all functions.\n\n\n\n\n\n","category":"function"},{"location":"man/tribe/#Tapestree.TRIBE.simulate_tribe","page":"TRIBE","title":"Tapestree.TRIBE.simulate_tribe","text":"simulate_tribe(X_initial::Float64,\n               nareas   ::Int64,\n               tree_file::String;\n               ωx       = 0.0,\n               σ²       = 0.5,\n               λ1       = 0.5,\n               λ0       = 0.2,\n               ω1       = 0.0,\n               ω0       = 0.0,\n               const_δt = 1e-4)\n\nSimulate tribe model.\n\n...\n\nArguments\n\nX_initial::Float64: trait starting value.\nnareas   ::Int64: number of areas.\ntree_file::String: full path to tree file.\nωx::Float64 = 0.0: simulated value of ω_x.\nσ²::Float64 = 0.5: simulated value of σ^2.\nω1::Float64 = 0.0: simulated value of ω_1.\nω0::Float64 = 0.0: simulated value of ω_0.\nλ1::Float64 = 0.5: simulated value of λ_1.\nλ0::Float64 = 0.2: simulated value of λ_0.\nconst_δt = 1e-4: # delta t used to approximate the simulation (lower values\n\nare more accurate but at a computation cost). ...\n\n\n\n\n\n","category":"function"},{"location":"#Tapestree.jl","page":"Home","title":"Tapestree.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Unravelling the evolutionary tapestry: Tapestree is a  Julia package of phylogenetic analyses of  diversification, trait and biogeographic dynamics.","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Tapestree currently holds the following phylogenetic models:","category":"page"},{"location":"","page":"Home","title":"Home","text":"\"Trait and Range Interspecific Biogeographic Evolution\" (TRIBE): which implements a joint model of trait evolution and biogeographic history as described in Quintero, I., & Landis, M. J. (2020). Interdependent phenotypic and biogeographic evolution driven by biotic interactions. Systematic biology, 69(4), 739-755. https://doi.org/10.1093/sysbio/syz082\n\"Environmental and State dependent Speciation and Extinction\" (ESSE): which implements a joint geographic and environmental model of diversification as described in Quintero, I., Landis, M. J., Jetz, W., & Morlon, H. (2022). The build-up of the present-day tropical diversity of tetrapods. Proceedings of the National Academy of Sciences, 2023. 120 (20) e2220672120. https://doi.org/10.1073/pnas.2220672120","category":"page"},{"location":"#Manual","page":"Home","title":"Manual","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n    \"man/installation.md\",\n    \"man/tribe.md\",\n    \"man/esse.md\",\n]\nDepth = 3","category":"page"},{"location":"#Contact","page":"Home","title":"Contact","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Contact Me","category":"page"},{"location":"man/installation/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"man/installation/#Installation-of-Julia","page":"Installation","title":"Installation of Julia","text":"","category":"section"},{"location":"man/installation/","page":"Installation","title":"Installation","text":"Download and follow the instructions here https://julialang.org/downloads/","category":"page"},{"location":"man/installation/#Installation-of-Tapestree","page":"Installation","title":"Installation of Tapestree","text":"","category":"section"},{"location":"man/installation/#Requirements:","page":"Installation","title":"Requirements:","text":"","category":"section"},{"location":"man/installation/","page":"Installation","title":"Installation","text":"Julia v1.x\nR installed\nR ape package installed.","category":"page"},{"location":"man/installation/","page":"Installation","title":"Installation","text":"Open julia and type the following","category":"page"},{"location":"man/installation/","page":"Installation","title":"Installation","text":"using Pkg\nPkg.add(\"Tapestree\")","category":"page"}]
}
