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
Tx = Union{sTxs, sTpe, sTfpe}




"""
    Union type for trait and rate data

Txs = Union{sTxs}
"""
Txs = Union{sTxs}



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

