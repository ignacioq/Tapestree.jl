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
              algorithm   ::String            = "flow",
              λpriors     ::Float64           = .1,
              μpriors     ::Float64           = .1,
              gpriors     ::Float64           = .1,
              lpriors     ::Float64           = .1,
              qpriors     ::Float64           = .1,
              βpriors     ::NTuple{2,Float64} = (0.0, 10.0),
              hpriors     ::Float64           = .1,
              optimal_w   ::Float64           = 0.8,
              screen_print::Int64             = 5,
              Eδt         ::Float64           = 0.01,
              ti          ::Float64           = 0.0,
              E0          ::Array{Float64,1}  = [0.0,0.0]) where {M,N}

  tv, ed, el, bts, x, y = 
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


  # make likelihood function
  if occursin(r"^[f|F][A-za-z]*", algorithm) 
    # prepare likelihood
    Gt, Et, X, p, fp, triads, lbts, ns, ned, nets, pupd, phid, nnps, nps,
    dcp, dcfp, pardic, k, h, ny, model, λevent!, rootll, assign_hidfacs! = 
      prepare_ll(cov_mod, tv, x, y, ed, el, bts, E0, h, 
        constraints = constraints, Eδt = Eδt, ti = ti)

    # make likelihood function
    llf = make_loglik(Gt, Et, X, triads, lbts, bts, ns, ned, nets, 
                      λevent!, rootll)
  
  elseif occursin(r"^[p|P][A-za-z]*", algorithm)




  end

  # create prior function
  lpf = make_lpf(pupd, phid, 
    λpriors, μpriors, gpriors, lpriors, qpriors, βpriors, hpriors, k, h, ny, model)

  # create posterior functions
  lhf = make_lhf(llf, lpf, assign_hidfacs!, dcp, dcfp, 
    Val(k), Val(h), Val(ny), Val(model))

  # run slice sampler
  R = slice_sampler(lhf, p, fp, nnps, nps, phid, length(pardic), 
                    niter, nthin, optimal_w, screen_print)

  write_ssr(R, pardic, out_file)

  return R
end






# A, b = ([0.008456227262140577 0.014485795428404051 0.012145153070847105 0.004708234087853492 0.008942365816778526 0.00904138405501543; 0.0042147950600335245 0.0073061273730702205 0.00606570920331544 0.0023751003834346605 0.004532853888917985 0.004572949142193423; 0.009289692221093623 0.015945932410316607 0.0133468426487302 0.005183158116424803 0.009852615190239045 0.00995792155757946; 0.005361271388844335 0.009191103055391867 0.007703789144818422 0.0030108936745224333 0.005722105780919311 0.005785810431822724; 0.0055876076148310945 0.009607174053385731 0.008033456003264123 0.0031510614422842276 0.005995864334005991 0.006059512601701913; 0.007087004688832629 0.012166041684552735 0.010186143879654186 0.003987653371043515 0.007582719719755497 0.007665318146827598], [0.005521766083171978, 0.002736066499272884, 0.1195389375401128, 0.15416990436058287, 0.12783101231802713, 0.5902023131988324])

# X0 = A\b
# X1 = inv(transpose(A)*A)*transpose(A)*b
# X2 = transpose(A)*inv(A*transpose(A))*b
# X3 = pinv(A)*b
# X4 = qr(A, Val(true))\b

# X0 ./ sum(X0)
# X1 ./ sum(X1)
# X2 ./ sum(X2)
# X3 ./ sum(X3)
# X4 ./ sum(X4)


# A0, b0 = ([0.025989634030602133 0.03980605523991394 0.03667344314301781 0.013072572116816318 0.023644912414139616 0.024467804993870426; 0.01152498644448025 0.02506472740537248 0.01716211441635794 0.0068890694083870545 0.014339365840497557 0.013809897779108583; 0.027109399142224887 0.043997653682059934 0.03855706266388582 0.014017561397518514 0.02598614874000256 0.026546549962941177; 0.013852854935276073 0.0223742406302683 0.019797386405038808 0.00809766204848449 0.015091408053782295 0.015460809764210345; 0.013692567202575662 0.024527581574083197 0.019889749814787402 0.008583877585081891 0.016622995986589863 0.016710369357708553; 0.0176387096155187 0.029846135049839323 0.025389830553426758 0.010642583723175551 0.0201881307414052 0.020502233924917165], [0.0002804778255828838, 0.00012411239624579215, 0.008438040780057738, 0.1331740022330538, 0.11216043799859805, 0.7458229287664618])

# A0 = rand(6,6)

# X0 = A0\b0
# X1 = inv(transpose(A0)*A0)*transpose(A0)*b0
# X2 = transpose(A0)*inv(A0*transpose(A0))*b0
# X3 = pinv(A0)*b0
# X4 = qr(A0, Val(true))\b0

# X0 ./ sum(X0)
# X1 ./ sum(X1)
# X2 ./ sum(X2)
# X3 ./ sum(X3)
# X4 ./ sum(X4)

# @benchmark transpose(A0)*inv(A0*transpose(A0))*b0
# @benchmark $A0\$b0
# @benchmark XX = qr($A0, Val(true))\$b0
# @benchmark ldiv!($XF, qr($A0, Val(true)), $b0)
# @benchmark XX = pinv(A0)*b0

# A0 = qr(A0, Val(true))

# @benchmark qr!(A0, Val(true))\b0



# A, b = ([0.14993439030644842 0.04796569175221402 0.14094180433496456 0.14828340968260248 0.011563658601371905 0.003921535102673917; 0.013918045884694807 0.23872901074895989 0.045336630927401325 0.021169786198734452 0.047046333675515355 0.001006605156550766; 0.09473178147720808 0.10616787583992532 0.1734281874085074 0.12103475834808254 0.022770388302737265 0.003914211425770966; 0.12561023452822626 0.03022835343063783 0.07742456056013272 0.16330665366673658 0.0095771244210081 0.0034440018538557855; 0.01116299375418916 0.07847204771288822 0.01599359600503739 0.01699768456996503 0.01782249492716639 0.0005086321857802653; 0.09838615384649697 0.04848029361769355 0.06438179644402427 0.1290901077740215 0.013066539086551998 0.0027781600237370826], [0.0007299442872090026, 0.00020257576229187223, 0.005509893342196102, 0.0167190070737138, 0.00684298444004459, 0.9699955950945446])
# X0 = zeros(6)
# Xp1 = zeros(6)
# ldiv!(X0, qr!(A, Val(true)), b)
# mul!(Xp1, A, X0)

# A = [0.015197449383890256 0.0032605392465147034 0.005562489148158521 0.015929783696091934 0.0008289996124442465 0.0001372771819319446; 0.015814211138033 0.003677930152752319 0.006060001951800452 0.016688651592304423 0.0009263125681241433 0.00014638267227162283; 0.021314801551459005 0.004805787989922332 0.008023464868278651 0.02243373260451748 0.0012146887805343292 0.0001954211225280912; 0.014505319251649964 0.0030844812432787214 0.005282814690529307 0.015193412250817216 0.0007850883007973472 0.0001306827834297119; 0.007115686546670057 0.001600735361302761 0.0026750671196537584 0.007487780140537985 0.0004047006638639422 6.519381510438745e-5; 0.01789913477057433 0.0038702842652995582 0.00657997787490324 0.018773502299578044 0.0009830980854043918 0.0001620537932287716]
# B = [0.039276293188995615 0.007740642767315528 0.013729126796817618 0.04089579805726122 0.001989201097427984 0.00034632545369448517; 0.03317832999609679 0.009836675095591722 0.014731524087025803 0.03584189763765932 0.00241690618506618 0.0003333300826230623; 0.043405319335169196 0.011202574864927204 0.017700147176751595 0.046252427724079476 0.0027899637385414423 0.0004156639147686485; 0.034903151938189075 0.006597042313086364 0.011916990492703785 0.036234912452986795 0.0017049115618393337 0.0003041594241250209; 0.010140336261288088 0.0026854884195277582 0.00419417711705097 0.010827868112034697 0.0006670181160248404 9.787702465644953e-5; 0.02580898913942561 0.005402019645383765 0.009310778559280444 0.026998734544822822 0.0013776748641390732 0.0002313908164451908]


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
  tree, bts = read_tree(tree_file, order = "postorder", branching_times = true)
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

  return tip_states, ed, el, bts, x, y
end





"""
  write_ssr(R       ::Array{Float64,2}, 
            pardic  ::Dict{String,Int64},
            out_file::String)

Write the samples from an MCMC data frame given a Dictionary of parameters.
"""
function write_ssr(R       ::Array{Float64,2}, 
                   pardic  ::Dict{String,Int64},
                   out_file::String)

  # column names
  col_nam = ["Iteration", "Posterior"]

  for (k,v) in sort!(collect(pardic), by = x -> x[2])
    push!(col_nam, k)
  end

  R = vcat(reshape(col_nam, 1, lastindex(col_nam)), R)

  writedlm(out_file*".log", R)
end





