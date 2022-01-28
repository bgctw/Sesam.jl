"""
    SystemParUpdater(d::Dict, sys::ModelingToolkit.AbstractSystem)

helps updating a subset of a parameter vector by another vector.     
The positions to update are specified during initialization by
providing the names or a dictionary.
The Callable can then be used to construct an updated Problem without
recreating the Problem from the System.

```julia
# given an ODESystem, sys, and associated ODEProblem, prob:

d_up = Dict(p3 => 0.3, p1 => 0.1) # a subset of the original parameters
upd! = SystemParUpdater(d_up,sys)
upd!(odeproblem_of_sys.p, d_up)
```
"""
struct SystemParUpdater{T}
    pos::T
end

function SystemParUpdater(d::Dict, sys::ModelingToolkit.AbstractSystem)
    SystemParUpdater(keys(d), sys)
end

function SystemParUpdater(d, sys::ModelingToolkit.AbstractSystem)
    pos = [findfirst(==(key), string.(parameters(sys))) for key in string.(d)]
    SystemParUpdater(pos)
end

function (upd::SystemParUpdater)(pall, dup::Dict)
    upd(pall, values(dup))
end
function (upd::SystemParUpdater)(pall, pup)
    pall[upd.pos] .= pup
    pall
end

# soft minimum - check if this reduces length of simplified formulas
# logaddexp does not work with Num
# smin(a,b) = -logaddexp(-a,-b)
# smin(args...) = -logaddexp(.-args)

# does not check x outside bounds to avoid comparison
logistic_unsafe(x;steepness=1e5) = 1 / (1 + exp(steepness * (-x)))
smin(a) = a
# smin(a,b) = logistic_unsafe(a-b) *b + logistic_unsafe(b-a) * a
smin(a,b) = (a>b) *b + (b>a) * a
# function smin(a,b,c) 
#     minab = smin(a,b)
#     smin(minab,c)
# end
# MAYBE implement with mapreduce for more arguments




