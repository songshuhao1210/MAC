function [] = myPlot(Plot_flag,freq,ri,m,N,L,EqS,unl,wnl,anl,cnl,r,path)

    %% Spectral coefficients to Physical space
    theta=linspace(0.0001,pi-0.0001,L+1);
    [Vel,Vor,Helicity]=Coef2Phy(ri,m,N,L,EqS,unl,wnl,r,theta);
    [MagB,CurrJ,HelB]=Coef2Phy(ri,m,N,L,1-EqS,anl,cnl,r,theta);
    %Mode=['m',num2str(m),'Fre',num2str(-imag(EigVal),'%.2e')]

    if Plot_flag == 0
        %%%% Equatorial plane
        phi=linspace(0,2*pi,L);
        Temp_eq_r=Vel(:,L/2+1,1);
        Temp_eq=real(Temp_eq_r)*cos(m*phi)-imag(Temp_eq_r)*sin(m*phi);
        xx=1.0*r*cos(phi);
        yy=1.0*r*sin(phi);
        figure
        hold on
        axis equal
        axis off
        pcolor(xx,yy,Temp_eq); shading interp;
        colormap('jet')
        xb=cos(linspace(0, 2*pi, 100));
        yb=sin(linspace(0, 2*pi, 100));
        plot(xb,yb,'k-','LineWidth',1.0);
        
        c = colorbar;
        ax = gca;
        axpos = ax.Position;
        c.Position(3) = 0.8*c.Position(3);
        ax.Position = axpos;
        c.Position(1) = 0.9;
        
        %fileformat = "D:/MAC/Le=%.1e_Ekman=%.1e_N2=%d/figure/omg=%.4f_%s.png";
        fileformat = strcat(path,"figure/omg=%.4f_%s.png");
        
        dd = freq;
        ss1 = 'EquatorialPlane';
        str = sprintf(fileformat,dd,ss1);
        pause(0.1)
        saveas(gcf,str)
        

        %%% CMB
        Data_r0=MagB(1,:,1);
        Data_CMB=real(Data_r0')*cos(m*phi)-imag(Data_r0')*sin(m*phi);
        lat=90-theta*180/pi;
        lon=phi*180/pi;
        figure('Position',[100 100 700 400],'color','w')
        axesm ('mollweid', 'Frame', 'on', 'Grid', 'on');
        %axesm ('hammer', 'Frame', 'on', 'Grid', 'on');
        pcolorm(lat,lon,Data_CMB);
        %colormap(mycolor);
        axis off
        hold on
        %h=colorbar;
        %title(['$\omega$=' num2str(-imag(EigVal),'%.2e')],'Interpreter','latex','FontSize',20)
        %export_fig([Mode,'CMB.png'],'-m 2')
        
        c = colorbar;
        ax = gca;
        axpos = ax.Position;
        c.Position(3) = 0.8*c.Position(3);
        ax.Position = axpos;
        c.Position(1) = 0.9;
        
        dd = freq;
        ss1 = 'CMB';
        str = sprintf(fileformat,dd,ss1);
        pause(0.1)
        saveas(gcf,str)


        %%% merid
        [tt,rr]=meshgrid(theta,r);
        x=rr.*sin(tt);
        y=rr.*cos(tt);

        figure('position',[100 100 400 400],'Color','w');
        subplot('Position',[0.15,0.12,0.35,0.75])
        axis off
        axis equal
        hold on
        pcolor(-x,y,real(Vel(:,:,3))); shading interp;
        
        c = colorbar('WestOutside');
        ax = gca;
        axpos = ax.Position;
        c.Position(3) = 0.8*c.Position(3);
        ax.Position = axpos;
        c.Position(1) = 0.1;
        
        ax1=gca;
        colormap(ax1,'jet');
        %plot([-ri, -1],[0,0],'k-');
        plot(-ri*sin(theta),ri*cos(theta),'k-');
        plot(-sin(theta),cos(theta),'k-');
        %xlim([-1,0])
        subplot('Position',[0.5,0.12,0.35,0.75])
        axis off
        axis equal
        hold on
        h2=pcolor(x,y,real(MagB(:,:,1))); shading interp;
        colorbar('EastOutside')
        ax2=gca;
        colormap(ax2);
        %colormap(ax2,mycolor);
        hold on
        plot([0, 0],[-ri,-1],'k-');
        plot([0, 0],[ri,1],'k-');
        plot(sin(theta),cos(theta),'k-');
        plot(ri*sin(theta),ri*cos(theta),'k-');
        
        %colorbar
        c = colorbar;
        ax = gca;
        axpos = ax.Position;
        c.Position(3) = 0.8*c.Position(3);
        ax.Position = axpos;
        c.Position(1) = 0.9;
        
        dd = freq;
        ss1 = 'merid';
        str = sprintf(fileformat,dd,ss1);
        pause(0.1)
        saveas(gcf,str)
        pause(0.1)

        %export_fig([Mode,'merid.png'],'-m 2')
        % figure('position',[100 100 300 600],'Color','w');
        % axis off
        % axis equal
        % hold on
        % pcolor(-x,y,real(Vel(:,:,3))); shading interp;
        % %colorbar('WestOutside')
        % ax1=gca;
        % colormap(mycolor);
        % export_fig Vel.png -transparent

        % figure('position',[100 100 300 600],'Color','w');
        % axis off
        % axis equal
        % hold on
        % pcolor(x,y,real(MagB(:,:,1))); shading interp;
        % %colorbar('WestOutside')
        % ax1=gca;
        % colormap(ax1,'jet');
        %export_fig Mag.png -transparent
    else
        subplot('position',[0.1,0.585,0.28,0.28])
        %%%% Equatorial plane
        phi=linspace(0,2*pi,L);
        Temp_eq_r=Vel(:,L/2+1,1);
        Temp_eq=real(Temp_eq_r)*cos(m*phi)-imag(Temp_eq_r)*sin(m*phi);
        xx=1.0*r*cos(phi);
        yy=1.0*r*sin(phi);
        hold on
        axis equal
        axis off
        pcolor(xx,yy,Temp_eq); shading interp;
        colormap('jet')
        xb=cos(linspace(0, 2*pi, 100));
        yb=sin(linspace(0, 2*pi, 100));
        plot(xb,yb,'k-','LineWidth',1.0);
        title("Equatorial plane")
        
        subplot(4,4,[9,10,11,12,13,14,15,16])
%       subplot('position',[0,0,0.75,0.5])
        %%% CMB
        Data_r0=MagB(1,:,1);
        Data_CMB=real(Data_r0')*cos(m*phi)-imag(Data_r0')*sin(m*phi);
        lat=90-theta*180/pi;
        lon=phi*180/pi;
%        figure('Position',[100 100 700 400],'color','w')
        axesm ('mollweid', 'Frame', 'on', 'Grid', 'on');
        %axesm ('hammer', 'Frame', 'on', 'Grid', 'on');
        pcolorm(lat,lon,Data_CMB);
        %colormap(mycolor);
        axis off
        hold on
        %h=colorbar;
        %title(['$\omega$=' num2str(-imag(EigVal),'%.2e')],'Interpreter','latex','FontSize',20)
        %export_fig([Mode,'CMB.png'],'-m 2')
        title("CMB")
       
        [tt,rr]=meshgrid(theta,r);
        x=rr.*sin(tt);
        y=rr.*cos(tt);

%       figure('position',[100 100 400 400],'Color','w');
%       subplot(4,4,[3,7])
        subplot('Position',[0.65,0.5,0.2,0.475])
        axis off
        axis equal
        hold on
        pcolor(-x,y,real(Vel(:,:,3))); shading interp;
        %colorbar('WestOutside')
        ax1=gca;
        colormap(ax1,'jet');
        %plot([-ri, -1],[0,0],'k-');
        plot(-ri*sin(theta),ri*cos(theta),'k-');
        plot(-sin(theta),cos(theta),'k-');
        %xlim([-1,0])
%        subplot(4,4,[4,8])
        subplot('Position',[0.65,0.5,0.2,0.475])
        axis off
        axis equal
        hold on
        h2=pcolor(x,y,real(MagB(:,:,1))); shading interp;
        %colorbar('EastOutside')
        ax2=gca;
        colormap(ax2);
        %colormap(ax2,mycolor);
        hold on
        plot([0, 0],[-ri,-1],'k-');
        plot([0, 0],[ri,1],'k-');
        plot(sin(theta),cos(theta),'k-');
        plot(ri*sin(theta),ri*cos(theta),'k-');
        title('merid')
        

    end
end
