macro autoinfiltrate(cond = true)
    pkgid = Base.PkgId(Base.UUID("5903a43b-9cc3-4c30-8d17-598619ec4e9b"), "Infiltrator")
    if !haskey(Base.loaded_modules, pkgid)
        try
            Base.eval(Main, :(using Infiltrator))
        catch err
            @error "Cannot load Infiltrator.jl. Make sure it is included in your environment stack."
        end
    end
    i = get(Base.loaded_modules, pkgid, nothing)
    lnn = LineNumberNode(__source__.line, __source__.file)
    if i === nothing
        return Expr(:macrocall, Symbol("@warn"), lnn, "Could not load Infiltrator.")
    end
    return Expr(:macrocall,
        Expr(:., i, QuoteNode(Symbol("@infiltrate"))), lnn, esc(cond))
end

function test_infiltrate()
    # only uncomment/use after precompilation, macro need to be defined in startup.jl
    #Main.@infiltrate_main 
    # works fine
    #if isdefined(Main, :Infiltrator); Main.infiltrate(@__MODULE__, Base.@locals, @__FILE__, @__LINE__); end
    # only works if @autoinfiltrate is defined in the package
    #@autoinfiltrate
    return "return from test_infiltrate()"
end
