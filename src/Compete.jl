
"""

Module for Compete.jl package

Ignacio Quintero 

t(-_-t)

22 June 2017

"""

module Compete

  export compete

  include("utils.jl")
  include("compete_utils.jl")
  include("cont_DA_prop.jl")
  include("data_initializer.jl")
  include("data_handling.jl")
  include("loglik_functions.jl")
  include("mcmc.jl")
  include("parameter_updates.jl")
  include("mcmc_burn.jl")
  include("proposal_functions.jl")
  include("compete_wrapper.jl")

end # module
