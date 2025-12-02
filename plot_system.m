function plot_system(xlist, tlist, Vlist, Uf_func, Vpred, string_params)
    close all
    figure
    n = size(Vlist,2)/2;
    p = plot(xlist,0,"-o", MarkerFaceColor="r",MarkerEdgeColor="r",MarkerSize=4, Color="k",LineWidth=2);
    if Vpred~=0
        hold on
        p1 = plot(xlist,0,"-o", MarkerFaceColor="c",MarkerEdgeColor="c",MarkerSize=4, Color="y",LineWidth=2);
    else
        l1 = xline(0);
    end
    axis([-0.5,max(xlist)+0.5,-0.5,0.5])
    tdiff = diff(tlist);
    txt = text(-0.3,0.3,string(tlist(1)));
    for j=1:length(tlist)-1
        tic()
        v = Vlist(j,1:n);
        if Vpred~=0
            r = Vpred(j,:);
            [~,ivmax] = max(abs(v));
            [~,irmax] = max(abs(r));
            set(p1,"xdata", xlist, "ydata",[0, r*v(irmax)/r(irmax), 0])
        else
            %short blurb showing how to find x-coord of tracking line
            %x = x-coord of tracking line, t = current time
            %c = wave speed, w = pulse width (in time), L = string length
            L = string_params.L;
            w = string_params.pulse_width;
            c = string_params.wave_speed;
            x = L-c*tlist(j)+.5*w*c;
            x = mod(x,2*L);
            if x > L
                x = 2*L - x;
            end
            set(l1,"value",x);
        end
        set(p,"xdata", xlist, "ydata",[0, v, Uf_func(tlist(j))])
        set(txt, 'string',string(tlist(j+1)))
        drawnow;
        slo=toc();
        pause(max(0,tdiff(j)-slo));
    end
end