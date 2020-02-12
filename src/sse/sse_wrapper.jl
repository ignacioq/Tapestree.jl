#=

Wrapper

Ignacio Quintero Mächler

t(-_-t)

September 26 2017

=#


"""
    ESSE(states_file ::String,
         tree_file   ::String,
         envdata_file::String,
         cov_mod     ::NTuple{M,String},
         out_file    ::String,
         h           ::Int64;
         constraints ::NTuple{N,String}  = (" ",),
         niter       ::Int64             = 10_000,
         nthin       ::Int64             = 10,
         λpriors     ::Float64           = .1,
         μpriors     ::Float64           = .1,
         gpriors     ::Float64           = .1,
         lpriors     ::Float64           = .1,
         qpriors     ::Float64           = .1,
         βpriors     ::NTuple{2,Float64} = (0.0, 10.0),
         hpriors     ::Float64           = .1,
         optimal_w   ::Float64           = 0.8,
         screen_print::Int64             = 5)

Wrapper for running a SSE model.
"""
function ESSE(states_file ::String,
              tree_file   ::String,
              envdata_file::String,
              cov_mod     ::NTuple{M,String},
              out_file    ::String,
              h           ::Int64;
              constraints ::NTuple{N,String}  = (" ",),
              niter       ::Int64             = 10_000,
              nthin       ::Int64             = 10,
              scale_y     ::NTuple{2,Bool}    = (true, false),
              λpriors     ::Float64           = .1,
              μpriors     ::Float64           = .1,
              gpriors     ::Float64           = .1,
              lpriors     ::Float64           = .1,
              qpriors     ::Float64           = .1,
              βpriors     ::NTuple{2,Float64} = (0.0, 10.0),
              hpriors     ::Float64           = .1,
              optimal_w   ::Float64           = 0.8,
              screen_print::Int64             = 5) where {M,N}

  tip_val, ed, el, x, y = 
    read_data_esse(states_file, tree_file, envdata_file)

  # if scale y
  if scale_y[1]
    # if scale each function separately or together
    if scale_y[2]
      ymin = minimum(y)
      ymax = maximum(y)
      for j in axes(y,2), i in axes(y,1)
        y[i,j] = (y[i,j] - ymin)/(ymax - ymin)
      end
    else
      ymin = minimum(y, dims = 1)
      ymax = maximum(y, dims = 1)
      for j in axes(y,2), i in axes(y,1)
        y[i,j] = (y[i,j] - ymin[j])/(ymax[j] - ymin[j])
      end
    end
  end

  R = slice_sampler(tip_val, ed, el, x, y, cov_mod, out_file, h,
        constraints  = constraints,
        niter        = niter,
        nthin        = nthin,
        λpriors      = λpriors,
        μpriors      = μpriors,
        gpriors      = gpriors,
        lpriors      = lpriors,
        qpriors      = qpriors,
        βpriors      = βpriors,
        hpriors      = hpriors,
        optimal_w    = optimal_w,
        screen_print = screen_print)

  return R
end





"""
    read_data_esse(states_file ::String, 
                   tree_file   ::String, 
                   envdata_file::String)

Process tree and state and environmental data file to run ESSE.
"""
function read_data_esse(states_file ::String, 
                        tree_file   ::String, 
                        envdata_file::String)

  # read tree in postorder and assign to objects
  tree = read_tree(tree_file, order = "postorder", branching_times = false)
  ntip = tree.nnod + 1
  ed   = tree.ed
  el   = tree.el
  tlab = tree.tlab

  # assign tip labels to edge numbers
  tip_labels = Dict{String,Integer}()
  ii = 0
  for i in Base.OneTo(size(ed,1))
    if ed[i,2] <= ntip
      ii += 1
      tip_labels[tlab[ii]] = ed[i,2]
    end
  end

  # read states text file
  data = DelimitedFiles.readdlm(states_file)

  if size(data,1) != ntip
    data = DelimitedFiles.readdlm(states_file, '\t', '\r')
  end

  if size(data,1) != ntip
    data = DelimitedFiles.readdlm(states_file, '\t', '\n')
  end

  if size(data,1) != ntip 
    error("Data file cannot be made of the right dimensions.\n Make sure the data file has the same number of rows as tips in the tree")
  end

  data_tlab    = convert(Array{String,1}, data[:,1])
  data_states  = convert(Array{Float64,2},  data[:,2:end])

  # create dictionary
  tip_states = Dict(tip_labels[val] => data_states[i,:] 
                   for (i,val) = enumerate(data_tlab))

  # process environmental data file
  envdata = DelimitedFiles.readdlm(envdata_file)

  x = envdata[:,1]
  y = envdata[:,2:end]

  return tip_states, ed, el, x, y
end






