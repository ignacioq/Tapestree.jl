#=

Ignacio Quintero Mächler

t(-_-t)

22 June 2017

=#

__precompile__()

"""

Tapestree package

"""

module Tapestree

  export compete
  export simulate_compete
  export compete_for_sims

  include("utils.jl")
  include("compete_utils.jl")
  include("cont_DA_prop.jl")
  include("disc_DA_prop.jl")
  include("data_initializer.jl")
  include("area_lineage_averages.jl")
  include("loglik_functions.jl")
  include("mcmc.jl")
  include("parameter_updates.jl")
  include("mcmc_burn.jl")
  include("proposal_functions.jl")
  include("compete_wrapper.jl")
  include("sim_utils.jl")

end # module Tapestree
