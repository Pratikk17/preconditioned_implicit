function Plot_Solution(Tout,Uout,flag,dim_prb)
Globals2D;
if flag==1
    Size_T = length(Tout);
    figure;    FigHandle = gcf;
    for k = 1:Size_T
        %Hx
        Hx = Uout(k,1:1/3*dim_prb); 
        Hx = reshape(Hx, Np, K);
        figure(1)
        PlotField2D(N,x,y,Hx);
        xlim([-1, 1]);        ylim([-1, 1]);        zlim([-1, 1]);
        grid on;         
        caxis manual;         caxis([-1, 1]);
        set(FigHandle, 'colormap', jet)
        set(FigHandle, 'Position', [200, 200, 800, 800]);
        set(FigHandle, 'PaperPositionMode', 'auto');
        title(['Time= ',num2str(Tout(k))]);
        colorbar('location', 'EastOutside'); 
        pause(0.01);

        %
        Hy = Uout(k,1/3*dim_prb+1:2/3*dim_prb); 
        Hy = reshape(Hy, Np, K);
        figure(2)
        PlotField2D(N,x,y,Hy);
        xlim([-1, 1]);        ylim([-1, 1]);        zlim([-1, 1]);
        grid on;         
        caxis manual;         caxis([-1, 1]);
        set(FigHandle, 'colormap', jet)
        set(FigHandle, 'Position', [200, 200, 800, 800]);
        set(FigHandle, 'PaperPositionMode', 'auto');
        title(['Time= ',num2str(Tout(k))]);
        colorbar('location', 'EastOutside'); 
        pause(0.01);

        Ez = Uout(k,2/3*dim_prb+1:end); 
        Ez = reshape(Ez, Np, K);
        figure(3)
        PlotField2D(N,x,y,Ez);
        xlim([-1, 1]);        ylim([-1, 1]);        zlim([-1, 1]);
        grid on;         
        caxis manual;         caxis([-1, 1]);
        set(FigHandle, 'colormap', jet)
        set(FigHandle, 'Position', [200, 200, 800, 800]);
        set(FigHandle, 'PaperPositionMode', 'auto');
        title(['Time= ',num2str(Tout(k))]);
        colorbar('location', 'EastOutside'); 
        pause(0.01);
    end
end
end

