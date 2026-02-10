#=

Union types

Ignacio Quintero Mächler

t(-_-t)

Created 19 01 2026
=#




"""
    Union type for label trees

Tlabel = Union{sT_label, sTf_label}
"""
Tlabel = Union{sT_label, sTf_label}




"""
    Union type for fossil data

iTf = Union{sTf_label, sTfbd, sTfbdX, iTfbd, iTfbdX}
"""
iTf = Union{sTf_label, sTfbd, cTfbd, acTfbd, iTfbd, sTfpe, sTxs}




"""
    Union type for unlabelled fossil data

uTf = Union{sTfbd, sTfbdX, iTfbd, iTfbdX}
"""
uTf = Union{sTfbd, cTfbd, acTfbd, iTfbd, sTxs}




"""
    Union type for gbm-bd data

iTbdU = Union{iTbd, iTfbd}
"""
cTbdU = Union{cTbd, cTfbd, acTfbd}




"""
    Union type for gbm-bd data

iTbdUX = Union{iTbdX, iTfbdX}
"""
iTbdU = Union{iTbd, iTfbd}




"""
    Union type for simple trait data

Tx = Union{sTpbx, sTbdx, sTfbdx}
"""
Tx = Union{sTxs, sTpe, sTfpe, iTxb}




"""
    Union type for trait and rate data

Txs = Union{sTxs}
"""
Txs = Union{sTxs, iTxb}



"""
    Union type for trait and rate data

Txs = Union{sTxs}
"""
Tpe = Union{sTpe, sTfpe}




"""
    Union type for asymmetrical trees

aT = Union{sTpe, sTfpe, acTfbd}
"""
aT = Union{sTpe, sTfpe, acTfbd}






"""
    sT_label(tree::T) where {T <: iTree}

Demotes a tree to `sT_label` and returns tip end traits.
"""
function sT_label(tree::T) where {T <: Tx}
  xs = Dict{String, Float64}()
  tr = _sT_label(tree, 0, xs)[1]
  return tr, xs
end


"""
    _sT_label(tree::T, i::Int64) where {T <: iTree}

Demotes a tree to `sT_label`.
"""
function _sT_label(tree::T, i::Int64, xs::Dict{String, Float64}) where {T <: Tx}
  if def1(tree)
    t1, i = _sT_label(tree.d1, i, xs)
    t2, i = _sT_label(tree.d2, i, xs)
    lab   = isdefined(tree, :l) ? label(tree) : ""
    tree  = sT_label(t1, t2, e(tree), lab)
  else
    i += 1
    lab = isdefined(tree, :l) ? label(tree) : string("t", i)
    xs[lab] = xv(tree)[end]
    tree = sT_label(e(tree), lab)
  end

  return tree, i
end

