#=

running tribe model from the command line

Ignacio Quintero Mächler

t(-_-t)

August 17 2017
=#



"""
    parse_commandline()

Parses command line arguments to be used in compete
"""
function parse_commandline()

  s = ArgParseSettings(prog        = "Tapestree",
                       description = "Inference of competition on trait evolution along biogeographic history",
                       version     = "0.1",
                       add_version = true)

  @add_arg_table s begin

    # optional arguments
    "--niter", "-i"
      help = "number of mcmc iterations"
      arg_type = Int64
      default = 500000

    "--nburn", "-b"
      help = "number of iterations in burning and tuning phase"
      arg_type = Int64
      default = 500000

    "--nthin", "-t"
      help = "sample every nthin iterations"
      arg_type = Int64
      default = 1000

    "--min_dt", "-d"
      help = "minimum percentage of discrete sampling"
      arg_type = Float64
      default = 0.01

    "--comp_x_prior"
      help = "gaussian prior for trait competition"
      arg_type = Tuple{Float64,Float64}
      default = 0.0, 10.0
      dest_name = "ωxprior"

    "--comp_gain_prior"
      help = "gaussian prior for effect of competition on area gain"
      arg_type = Tuple{Float64,Float64}
      default = 0.0, 10.0
      dest_name = "ω1prior"

    "--comp_loss_prior"
      help = "gaussian prior for effect of competition on area loss"
      arg_type = Tuple{Float64,Float64}
      default = 0.0, 10.0
      dest_name = "ω0prior"

    "--sigma_prior"
      help = "exponential prior for rate of evolution"
      arg_type = Float64
      default = 0.1
      dest_name = "σ²prior"

    "--lambda_prior"
      help = "exponential prior for gain and loss rates"
      arg_type = Float64
      default = 0.1
      dest_name = "λprior"

    "--stbrl"
      help = "stem branch length"
      arg_type = Float64
      default = 1.0

    # flags
    "--fix_biocomp", "-f"
      help = "fix to 0 competition effect on gain and losses"
      action = :store_true
      dest_name = "fix_ω1_ω0"

    # required arguments 
    "tree_file"
      help = "full path to phylogenetic tree file"
      required = true
    "data_file"
      help = "full path to data file"
      required = true
    "out_file"
      help = "full path to write file"
      required = true

  end

  parsed_args = parse_args(s, as_symbols = true)

  println("Parsed args:")
  for (arg,val) in parsed_args
      println("  $arg  =>  $val")
  end

  return parsed_args
end

# parse commands
comp_args = parse_commandline()

using Tapestree

# run compete
tribe(comp_args[:tree_file], 
      comp_args[:data_file], 
      comp_args[:out_file];
      filter((u,v) -> u != :tree_file && 
                      u != :data_file && 
                      u != :out_file, comp_args)...)


