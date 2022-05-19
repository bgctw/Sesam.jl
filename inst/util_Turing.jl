function Chains_section(a::AbstractArray, section::Symbol)
    varnames = [string(section)*"["*string(i)*"]" for i in 1:size(a,2)]
    Chains_section(a, section, varnames)
end
function Chains_section(a::AbstractArray, section::Symbol, varnames)
    Chains(a, varnames, Dict(section => varnames))
end

"""
    add_vars(fgetvars, chn, section::Symbol, var_ret, var_arg = names(chn,[:parameters]))

Extend Chain `chn` by section consisting of variables var_ret. The variables are
computed by function fgetvars(args), where args is a subset of a single sample in chain.
The subset is taken by extracting variables var_arg.

TODO example, testcase, used in exploreP.add_obs
"""
function add_vars(fgetvars::Function, chn, section::Symbol, var_ret, var_arg = names(chn,[:parameters]))
    # var_ret = symbol.(var_obs)
    # var_arg = names(chn,[:parameters, :u0])
    #tmp_vars = mapslices(identity, chn[:,index_vars,:].value, dims=(2))
    tmp_vars = mapslices(fgetvars, chn[:,var_arg,:].value, dims=(2))
    chn_obs = Chains(tmp_vars, symbol.(var_ret), Dict(section => var_ret))
    if chains(chn) != chains(chn_obs)
        chn = Chains(chn.value, names(chn), chn.name_map)
    end
    chn_ext = hcat(chn, chn_obs)
end

function chainscat_resample(chns) 
    # reduce to subset of variables
    par_names = names.(chns)
    par_names_int = reduce(intersect, par_names)
    i_chn_morevars = length.(par_names) .> length(par_names_int)
    chns2 = map(chns, i_chn_morevars) do chn, is_more
        is_more ? chn[:,par_names_int,:] : chn
    end
    # reduce to same minimum sample number
    fsubsample = (chn) -> Chains(chn[rand(1:size(chn,1), n_min),:,:].value, names(chn), chn.name_map)
    n_min, n_max = extrema(map(x -> size(x,1), chns))
    chns3 = map(chns2) do chn
        size(chn,1) == n_min ? chn : fsubsample(chn)
    end
    chn_exti = chainscat(chns3...)
end

"""

Construct chain of given parameters of original chain. This helps to 
plot variables that are not in the :parameters section
"""
chains_par(chns, par_names) = Chains(chns[:,par_names,:].value, par_names)
    
