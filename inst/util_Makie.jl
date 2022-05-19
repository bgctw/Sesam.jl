import CairoMakie as CM
import KernelDensity

# import Makie.MakieLayout: defined_interval, defaultlimits
# import Makie: inverse_transform    

# function traceplot_params(chns, pars=names(chns, :parameters); fig = CM.Figure(; resolution=(1_000, 800)), column = 1)
#     n_chains = length(chains(chns))
#     n_samples = length(chns)
#     #i,param = first(enumerate(pars))
#     for (i, param) in enumerate(pars)
#         ax = CM.Axis(fig[i, column]; ylabel=string(param))
#         # i_chain = 1
#         for i_chain in 1:n_chains
#             _values = chns[:, param, i_chain]
#             CM.lines!(ax, 1:n_samples, _values; label=string(i_chain))
#         end
#         CM.hideydecorations!(ax; label=false)
#         if i < length(pars)
#             CM.hidexdecorations!(ax; grid=false)
#         else
#             ax.xlabel = "Iteration"
#         end
#     end
#     fig    
# end

function density_params(chns, pars=names(chns, :parameters), cfg::TWP.MakieConfig=TWP.MakieConfig(); fig = pdf_figure(cm2inch.((8.3,8.3/1.618)), MakieConfig()), column = 1, xlims=nothing, labels=nothing, colors = nothing, ylabels = nothing, normalize = false, kwargs_axis = repeat([()],length(pars)), kwargs...)
    n_chains = size(chns,3)
    n_samples = length(chns)
    labels_ch = isnothing(labels) ? string.(1:n_chains) : labels
    ylabels = isnothing(ylabels) ? string.(pars) : ylabels
    !isnothing(xlims) && (length(xlims) != length(pars)) && error(
        "Expected length(xlims)=$(length(xlims)) (each a Tuple or nothing) to be length(pars)=$(length(pars))")
    for (i, param) in enumerate(pars)
        ax = CM.Axis(fig[i, column]; ylabel=ylabels[i], kwargs_axis[i]...)
        if isnothing(colors)
            colors = ax.palette.color[]
        end
        for i_chain in 1:n_chains
            _values = chns[:, param, i_chain]
            col = colors[i_chain]
            if normalize
                k = KernelDensity.kde(_values)
                md = maximum(k.density)
                CM.lines!(ax, k.x, k.density ./ md; label=labels_ch[i_chain], color = col, kwargs...)
            else
                CM.density!(ax, _values; label=labels_ch[i_chain], color = (col, 0.3), strokecolor = col, strokewidth = 1, 
                #strokearound = true,
                kwargs...)
            end
        end
        xlim = passnothing(getindex)(xlims, i)
        !isnothing(xlim) && CM.xlims!(ax, xlim)
    #CM.hideydecorations!(ax,  ticklabels=false, ticks=false, grid=false)
        CM.hideydecorations!(ax, label=false, ticklabels=true)
        # if i < length(params)
        #     hidexdecorations!(ax; grid=false)
        # else
        #     ax.xlabel = "Parameter estimate"
        # end
    end
    # axes = [only(contents(fig[i, 2])) for i in 1:length(params)]
    # linkxaxes!(axes...)
    #axislegend(only(CM.contents(fig[2, column])))
    fig    
end


i_tmp = () -> begin
    chns = chn_exti;
    pars=names(chns, :parameters)
    pars=[:pl₊β_Ni0,:pl₊β_Pi0]
    column = 1
    #fig = CM.Figure(; resolution=(1_000, 800));
    fig,ax = pdf_figure_axis(cm2inch.((8.3,7.0)));  delete!(fig.current_axis[])
    density_params(chn_exti, [:s₊a_E, :s₊i_BP]; fig, xlims=[(0,200),(0,30)], labels=scen_df.site)
    density_params(chn_exti, [:pl₊β_Ni0,:pl₊β_Pi0]; fig, labels=scen_df.site)
    CM.axislegend(only(CM.contents(fig[1, column])))
    fig
end

# struct UnitMultiplier{T}
#     factor::T
# end
# function (m::UnitMultiplier)(x)
#     m.factor * x
# end
# Makie.MakieLayout.defined_interval(::UnitMultiplier) = Makie.OpenInterval(-Inf, Inf)
# Makie.MakieLayout.defaultlimits(::UnitMultiplier) = Makie.MakieLayout.defaultlimits(identity)
# Makie.inverse_transform(m::UnitMultiplier) = (x) -> x/m.factor

