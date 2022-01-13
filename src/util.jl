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



