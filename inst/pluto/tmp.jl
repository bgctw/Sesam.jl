### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 884883c1-c2bf-4209-b783-2e57428e4389
if occursin("Sesam", Base.current_project())
	import Pkg
	# activate the shared project environment
	Pkg.project().path != Base.current_project() && 
		Pkg.activate(Base.current_project())
	Pkg.develop("MTKHelpers")
	# instantiate, i.e. make sure that all packages are downloaded
	Pkg.instantiate()
	using PlutoLinks: @revise
	using Sesam, MTKHelpers
	#@revise using Sesam
	#@revise using MTKHelpers
else
	using Sesam, MTKHelpers
end

# ╔═╡ 27f3db92-b75d-11ec-00ad-3bcc21f6f62d
Base.current_project()


# ╔═╡ cb96aaa0-dd01-42c0-816b-cd3a32c9966c
pwd()

# ╔═╡ e15fcb5c-f5d2-4ce8-9d65-5fe0f0959cab
# begin
# 	import Pkg
# 	Pkg.activate(Base.current_project())
# 	Pkg.project().path == Base.current_project()
# end

# ╔═╡ b2336889-8952-4e52-8aa2-d95efcb78dc9
Pkg.project().path

# ╔═╡ Cell order:
# ╠═27f3db92-b75d-11ec-00ad-3bcc21f6f62d
# ╠═cb96aaa0-dd01-42c0-816b-cd3a32c9966c
# ╠═e15fcb5c-f5d2-4ce8-9d65-5fe0f0959cab
# ╠═884883c1-c2bf-4209-b783-2e57428e4389
# ╠═b2336889-8952-4e52-8aa2-d95efcb78dc9
