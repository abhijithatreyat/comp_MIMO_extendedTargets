function plot_NOMP_results(PlotResultsFlag, ffreq_vect, doppler_vect, ...
    ffreq_to_range, doppler_to_speed, RangeList, DopplerList, gainListRD,Rmin,Rmax)
if PlotResultsFlag
    fig_NOMP = figure();
    xlabel('range (m)'); ylabel('radial speed (m/s)');
    title('true (x) and estimated (o) targets');
    
    for fig_i = 1:2*2
        figure(fig_NOMP);
        subplot(2,2,fig_i);
        xlabel('range (m)'); ylabel('radial speed (m/s)');
        title(['Beacon',num2str(fig_i)]);
        grid on;
        plot(ffreq_vect*ffreq_to_range,doppler_vect*doppler_to_speed,'rx','LineWidth',.75);
        hold on;
        xlim([Rmin-1,Rmax+1]); ylim([-1,1]*pi*doppler_to_speed) % Xlim 

        for val = 1:size(RangeList,1)
            plot(RangeList(val,fig_i)',DopplerList(val,fig_i)','bo','LineWidth',.75,'MarkerSize', ...
                1+floor(abs(gainListRD(val,fig_i)')) );
            hold on;
        end
        hold off;
    end
end
end