function color_legend_texts(leg)
    """Color legend texts based on color of corresponding lines"""
    for (line, txt) in zip(leg.get_lines(), leg.get_texts())
        txt.set_color(line.get_color())
    end
end

function color_list(time_run)
    cmap = get_cmap("rainbow")
    color_array = cmap.(LinRange(0, 1, time_run))
    return color_array
end
