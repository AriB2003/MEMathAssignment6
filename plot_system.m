%animate and plot the motion of a vibrating string over time.
%INPUTS
%xlist: vector of x-positions
%tlist: vector of time values where state is sampled
%Vlist: matrix of system states; each row is [U; dU/dt] at a given time
%       where U is the displacement of the n internal masses
%Uf_func: function handle giving displacement of the right boundary at time t
%Vpred: predicted mode shape (optional). If zero, a tracking line is shown
%string_params: struct containing system parameters for reference/auxiliary logic
%   string_params.L: string length
%   string_params.wave_speed: wave propagation speed
%   string_params.pulse_width: temporal pulse width (used when Vpred=0)
%tit: title string for the animation figure
%writerObj: VideoWriter object used to record frames to a video file
function plot_system(xlist, tlist, Vlist, Uf_func, Vpred, string_params, tit, writerObj)
    close all
    fig1 = figure();
    n = size(Vlist,2)/2; % number of masses (state vector = [U; dU/dt])

    % initialize main plot
    p = plot(xlist,zeros(size(xlist)),"-o", MarkerFaceColor="r",MarkerEdgeColor="r",MarkerSize=4, Color="k",LineWidth=2, DisplayName="Modeled");
    
    % optional tracking line
    if Vpred~=0
        hold on
        p1 = plot(xlist,zeros(size(xlist)),"-o", MarkerFaceColor="c",MarkerEdgeColor="c",MarkerSize=4, Color="y",LineWidth=2, DisplayName="Predicted");
    else
        l1 = xline(0, DisplayName="Wavespeed Tracker");
    end


    legend()
    title(tit)
    axis([-0.5,max(xlist)+0.5,-0.5,0.5])
    tdiff = diff(tlist); % time step durations for animation pacing
    txt = text(-0.3,0.3,string(tlist(1)));  % timestamp label
    
    % main animation loop
    for j=1:length(tlist)-1
        tic()
        v = Vlist(j+1,1:n);
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
        set(p,"xdata", xlist, "ydata",[0, v, Uf_func(tlist(j+1))])
        set(txt, 'string',string(tlist(j+1)))
        drawnow;
        %capture a frame (what is currently plotted)
        current_frame = getframe(fig1);
        %write the frame to the video
        writeVideo(writerObj,current_frame);
        slo=toc();
        pause(max(0,tdiff(j)-slo));
    end
    for i = 1:60
        writeVideo(writerObj,current_frame);
    end
end