function plotError(t :: Array{Float64,1},
                    err :: Array{Float64,1},
                    error_norm :: String;
                    error_str :: String = "relative error ",
                    plotcolor :: String = "blue",
                    style :: String = "-",
                    labelstr :: String = "",
                    file_name :: String = "")
    
    font_1 = Dict("family"=>"serif",
                  "color"=>"darkred",
                  "weight"=>"normal",
                  "size"=>14)

    # create a figure and save its handle

    fig = gcf()
    fig[:set_size_inches](15,8, true)

    ax = subplot(111)
    ax[:set_ylim]([0;2.0])
    ax[:spines]["top"][:set_color]("none")
    setp(ax[:get_yticklabels](), fontsize = 14, color = "darkred") 


    pl_1 = ax[:plot](t, err)
    setp(pl_1, linestyle = style,
         marker = "None",
         linewidth = 2,
         color = plotcolor,
         label = labelstr) 


    legend(loc = "upper right", fancybox = "true") 
    ylabel(error_str, fontdict = font_1)
    xlabel("time", fontdict = font_1)

    if file_name != ""
        matplotlib[:pyplot][:savefig](string("data/", file_name, "norm_$(error_norm)",".pdf"), transparent=true, format="pdf", dpi=1000)
        close()
    end 

end