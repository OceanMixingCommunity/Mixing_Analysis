% plot the Kt over Kr for andy

%  Kt/krho (using eqs 1 +2)
% KT =1/2 * chi / < dT=dz >^2 (1)
% Kr =0.2 * eps / N^2 (2)
colz = cbrewer('qual','Paired',26);

for i=[1:3,5,10]
    if i==3;

        tnm = 'Natre';
        load('../../data/natre/natre_hrp_for_map.mat')
        
        K = 0.2 * hrp.eps ./ hrp.N2; 
        presK = hrp.pres;
        
        KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
    
    elseif i ==1;

        tnm = 'BBTRE (smooth)';
        load('../../data/bbtre/bbtre_hrp_for_map.mat')
        hrp96=hrp;  %[1:20,71:75]
        load('../../data/bbtre/bbtre97_hrp_for_map.mat')
        hrp97=hrp;  % 1:36    
        
        hrp_bbtre_smooth = 0.2*([hrp96.eps(:, [2:5,7:20,71:74])])./([hrp96.N2(:, [2:5,7:20,71:74])]);
        
        K = hrp_bbtre_smooth; 
        presK=hrp96.pres(:, [2:5,7:20,71:74]);
        
        KT = 0.5 * hrp96.chit(1:end-1,[2:5,7:20,71:74]) ./ (diff(hrp96.temp(:,[2:5,7:20,71:74]))./diff(hrp96.pres(:,[2:5,7:20,71:74]))).^2;
    
    elseif i ==2;

        tnm = 'BBTRE (rough)';
        load('../../data/bbtre/bbtre_hrp_for_map.mat')
        hrp96=hrp;
        load('../../data/bbtre/bbtre97_hrp_for_map.mat')
        hrp97=hrp;
        
        hrp_bbtre_rough = 0.2*([hrp96.eps(1:10815,[1,6,21:70,75]) hrp97.eps(:,1:end)])./([hrp96.N2(1:10815,[1,6,21:70,75]) hrp97.N2(:,1:end)]);
        K = hrp_bbtre_rough; 
        presK = [hrp96.pres(1:10815,[1,6,21:70,75]) hrp97.pres(:,1:end)];
            
        KT = 0.5 * ([hrp96.chit(1:10815-1,[1,6,21:70,75]) hrp97.chit(1:end-1,1:end)])./  (diff([hrp96.temp(1:10815,[1,6,21:70,75]) hrp97.temp(:,1:end)])./diff([hrp96.pres(1:10815,[1,6,21:70,75]) hrp97.pres(:,1:end)])).^2;

    elseif i ==4;

        tnm = 'Ladder';
        load('../../data/ladder/ladder_hrp_for_map.mat')
        K = 0.2 * hrp.eps ./ hrp.N2; 
        presK = hrp.pres;
        
        KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
    elseif i ==5;
        tnm = 'Graviluck';

        % only plot down to 2216m since below that one profile exists
        load('../../data/graviluck/graviluck_hrp_for_map.mat')
        K = 0.2 * hrp.eps ./ hrp.N2; 
        presK = hrp.pres;
        
        KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
    elseif i ==6;

        tnm= 'Fieberling';
        load('../../data/fieberling/fieberling_hrp_for_map.mat')
        K = 0.2 * hrp.eps ./ hrp.N2; 
        presK = hrp.pres;
        
        KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
    elseif i ==7;

        tnm= 'Dimes_{DP}';
        load('../../data/dimes/dimes_hrp_for_map.mat')
        K = 0.2 * hrp.eps ./ hrp.N2; 
        presK = hrp.pres;
        
        KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
    elseif i ==8;

        tnm= 'Dimes_{West}';
        load('../../data/dimes/dimes_pacific_hrp_for_map.mat')
        K = 0.2 * hrp.eps ./ hrp.N2; 
        presK = hrp.pres;
        
        KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
    elseif i ==9;

        tnm= 'Toto';
        load('../../data/toto/toto_hrp_for_map.mat')
        K = 0.2 * hrp.eps ./ hrp.N2; 
        presK = hrp.pres;
        
        KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
    elseif i ==10;
        tnm = 'Geotraces';
        load('../../data/geotraces/geotraces_hrp_for_map.mat');
        K = 0.2 * hrp.eps ./ hrp.N2; 
        presK = hrp.pres;
        
        KT = 0.5 *  hrp.chi(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
        
    end
    
    figure('paperposition',[0 0 11 4],'color','w'); wysiwyg; 
    ax=subplot(221);
    pcolor(1:size(presK,2),nanmean(presK,2)',log10(abs(K)));
    axis ij; colorbar; 
    caxis([-8 -4])
    ylabel('Depth (m)');
    xlabel('log_{10} K_{\rho} (m^2 s^{-1})')
    grid on; shading flat; 
    title(tnm)
    colormap(ax,cbrewer('seq','Blues',11))

    ax=subplot(222);
    pcolor(1:size(presK,2),nanmean(presK(1:end-1,:),2)',log10(KT));
    axis ij; colorbar; 
    set(gca,'yticklabel',[])
    caxis([-8 -4])
    xlabel('log_{10} K_{T} (m^2 s^{-1})')
    grid on; shading flat; 
    colormap(ax,cbrewer('seq','Blues',11))

    ax2=subplot(223);
    pcolor(1:size(presK,2),nanmean(presK(1:end-1,:),2)',log10(KT./abs(K(1:end-1,:))));
    axis ij; colorbar; 
    caxis([-2 2])
    ylabel('Depth (m)')
    %set(gca,'yticklabel',[])
    xlabel('log_{10} K_{T} / K_{\rho} ')
    grid on; shading flat; 
    colormap(ax2,parula)

    fullratio = log10(KT./abs(K(1:end-1,:)));
    fullratio = fullratio(:);
    edges = [-5:.2:5];
    [n]=histcounts(fullratio,edges);
    
    ax2=subplot(224);
    bar(edges(1:end-1)+diff(edges(1:2))/2,n);
    xlabel('log10(K_T/K_{\rho})')
    ylabel('number of counts');
    set(gca,'xtick',[-5:1:5])
    grid on; 
    
    export_fig('-dpng','-r100',['KToverK_' tnm '.png']);

end
