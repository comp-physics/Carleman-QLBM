using PyPlot

rc("font",family="sans serif")
# Disable LaTeX rendering for CI compatibility - use Unicode characters instead
rc("text",usetex="False")
rc("font", size=16)

function fm_plot(fT, n_time, nm_sym, nm_color, nm_label)
    time = range_step(1, n_time, 1)

    row = 1
    col = size(fT)[1]

   figure()
#    subplots_adjust(left=0.06, right=.95, top=.95, bottom=0.07, hspace = 0.2, wspace=0.22)
    for m=1:col
        subplot(row, col, m)
        plot(time, fT[m, :], marker = nm_sym, markersize = 4, color = nm_color, label = nm_label)
        ylabel("f_$m")
        xlabel("Time")
    end
end

function plot_CLBM_LBM(fT, VT_f, n_time, nm_color, nm_label, nm_legend)
    time = range_step(1, n_time, 1)

#    close("all")
#    figure(figsize=(10, 8))
    row = 1
    col = size(fT)[1]
#    subplots_adjust(left=0.06, right=.95, top=.95, bottom=0.07, hspace = 0.2, wspace=0.22)
    for m=1:col
        subplot(row, col, m)
        plot(time, fT[m, :], "ok", label = nm_legend)
        plot(time, VT_f[m, :], "+",  color = nm_color,  label = nm_label)
        ylabel("f_$m")
        xlabel("Time")
    end
    legend(loc= "best")
end

function plot_CLBM_LBM_diff(fT, VT_f, n_time, nm_color, nm_label, nm_legend)
    time = range_step(1, n_time, 1)

    relative_error = abs.(VT_f .- fT) ./ fT 
    y_lim = [1.e-16, 2.e-14]
#    row = 1
#    col = size(fT)[1]
#    for m=1:col
    row = size(fT)[1]
    col = 1
    for m = 1:row
        subplot(row, col, m)
        semilogy(time, relative_error[m, :], linestyle = "-", marker = "o", markersize = 4, color = nm_color,  label = nm_label)
       # ylabel(string("Relative error of f", m))
       ylabel("ε_CLBM_$m")
#        xlabel("Time step")
        ylim(y_lim)
    end
#    legend(loc= "best")
end
