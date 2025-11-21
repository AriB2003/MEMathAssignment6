function plot_system(xlist, tlist, Vlist, Uf_func)
    close all
    figure
    n = size(Vlist,2)/2;
    p = plot(xlist,0,"-o", MarkerFaceColor="r",MarkerEdgeColor="r",MarkerSize=4, Color="k",LineWidth=2);
    axis([-0.5,max(xlist)+0.5,-0.5,0.5])
    tdiff = diff(tlist);
    for j=1:length(tlist)-1
        set(p,"xdata", xlist, "ydata",[0, Vlist(j,1:n), Uf_func(tlist(j))])
        pause(tdiff(j));
    end
end