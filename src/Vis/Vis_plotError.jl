function plotErrorL2(error_std :: Array{Float64}, error_ms :: Array{Float64},
						t_max :: Float64, test_name = [])

	t = linspace(0, t_max, length(error_std))

	font_1 = Dict("family"=>"serif",
                  "color"=>"black",
                  "weight"=>"normal",
                  "size"=>14)    
    
    fig = figure("L2 error", figsize = (15,12)) 

    ax_1 = subplot(111)
    ax_1[:spines]["top"][:set_color]("none")
    ax_1[:spines]["right"][:set_color]("none")

    # ax_1[:set_ylim](y_lim)
    
    ax_1[:yaxis][:set_ticks_position]("left")
    ax_1[:xaxis][:set_ticks_position]("bottom")
    
    setp(ax_1[:get_yticklabels](), fontsize = 14, color = "black")
    setp(ax_1[:get_xticklabels](), fontsize = 14, color = "black")
    #setp(ax_1[:ticklabel_format](axis = "y", style = "sci", scilimits = [-2, 2]))

    pl_1 = ax_1[:plot](t, error_std)
    setp(pl_1, linestyle = "-",
         marker = "None",
         linewidth = 2, color = "red",
         label = "standard FEM (low res)")

    pl_2 = ax_1[:plot](t, error_ms)
    setp(pl_2, linestyle = "-",
         marker = "None",
         linewidth = 2, color = "blue",
         label = "ms-reconstruction")
    
    legend(loc = "upper right", fontsize = 14, fancybox = "true") 

    title("L2-error ")
    xlabel("time", fontdict = font_1)
    ylabel("error", fontdict = font_1)

    if test_name != []
        matplotlib[:pyplot][:savefig](string(test_name, "---L2-error ", ".pdf"), format="pdf", dpi=1000)
        close()
    end
    
    return nothing
end