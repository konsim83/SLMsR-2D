ENV["JULIA_NUM_THREADS"] = 4
const PLOTS_DEFAULTS = Dict(:theme => :wong)

import Pkg

let
	pkgs = ["Revise", 
			"OhMyREPL",
			"Rebugger",
			"Statistics",
			"SpecialFunctions",
			"PyPlot",
			"PyCall",
			"WriteVTK",
			"ProgressMeter",
			"NearestNeighbors",
			"BinDeps",
			"DifferentialEquations",
			"TriangleMesh"]
	for pkg in pkgs
		if Base.find_package(pkg) === nothing
		    Pkg.add(pkg)
		end
	end
end

try
	itr = walkdir(pwd())
	(root, dirs, files) = first(itr)

	push!(LOAD_PATH, root)

	for dir in dirs
		if ~(dir[1]=='.')
		  push!(LOAD_PATH, string(root, "/", dir))
		end	  

		# If we have a src and/or test directory add also this.
		if (dir == "src" || dir == "test" )
			itr_sub = walkdir(string(root, "/", dir))
			(root_sub, dirs_sub, files_sub) = first(itr_sub)
			for dir_sub in dirs_sub
				push!(LOAD_PATH, string(root_sub, "/", dir_sub))
			end
		end
	end	
catch
	@warn "Problem adding subdirectories to LOAD_PATH."
end

try
    @eval using Revise
    # Turn on Revise's automatic-evaluation behavior
    Revise.async_steal_repl_backend()

    @eval using Rebugger
    # Activate Rebugger's key bindings
    atreplinit(Rebugger.repl_init)

    using OhMyREPL
    colorscheme!("Monokai24bit")
catch err
    @warn "Could not load startup packages."
end

println("--------------------------------------------")
@info "The following paths were added to LOAD_PATH:"
display(LOAD_PATH)
println("\n--------------------------------------------")
