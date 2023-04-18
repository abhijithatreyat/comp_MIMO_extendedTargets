function [target_bins, detected_range] = extract_targets(N_target,power_bins, extra_bins,...
    power_residue, N_symb, oversampling_symb, N_chirp, oversampling_chirp, ...
    y_residue, y_rd_residue, N_beacons, ffreq_to_range, ts, PlotResultsFlag,...
    fig_0, fig_2, num_subplot_rows, num_subplot_cols, Rmin, Rmax, doppler_to_speed,...
    range_axis, speed_axis, min_val, peak_power)

% EXTRACT_TARGETS identifies targets from a range-doppler signal power residue
detected_range = [];

for i_target = 1:N_target+extra_bins
    [max_v,ind_row_vect] = max(power_residue);
    [~,ind_col] = max(max_v);
    ind_row = ind_row_vect(ind_col);
    target_bins(i_target,:) = [ind_row,ind_col]; % row and column index of 
    % largest peak in range-doppler signal
    
    w_range = 2*pi*(ind_col-1)/N_symb/oversampling_symb; % equivalent range frequency 
    % of largest peak
    w_doppler = angle(exp(-1i*2*pi*(ind_row-1)/N_chirp/oversampling_chirp)); % equivalent
    % doppler frequency of largest peak (identified target)
    y0 = reshape(exp(-1i*w_doppler*(0:N_chirp-1)')* ...
        exp(1i*w_range*(0:N_symb-1)),1,N_chirp,N_symb); % time domain response of 
                                                        %   identified target
    for i_beacon = 1:N_beacons % extract the identified target's response 
        % from signal residue in fast-slow time and range-doppler domain
        y_residue(i_beacon,:,:) = y_residue(i_beacon,:,:) - ...
            y_rd_residue(i_beacon,ind_row,ind_col)/N_symb/N_chirp*y0;
        y_rd_residue(i_beacon,:,:) = reshape(fft2(reshape(y_residue(i_beacon,:,:), ...
            N_chirp,N_symb),N_chirp*oversampling_chirp,N_symb*oversampling_symb),1, ...
            N_chirp*oversampling_chirp,N_symb*oversampling_symb);
    end
    power_residue = reshape(sum(abs(y_rd_residue).^2,1), ...
        N_chirp*oversampling_chirp,N_symb*oversampling_symb); % update aggregate power
    % residue
    detected_range = [ detected_range ,(ind_col-1)*2*pi/N_symb/oversampling_symb*ffreq_to_range/ts];
    if PlotResultsFlag
        if i_target <= num_subplot_rows*num_subplot_cols
            figure(fig_2); subplot(num_subplot_rows,num_subplot_cols,i_target);
            h = surf(range_axis,speed_axis,10*log10(power_residue+min_val));
            set(h,'edgecolor','none'); view(2);
            % xlim([0,N_symb*oversampling_symb]); ylim([0,N_chirp*oversampling_chirp]);
            xlim([Rmin,Rmax]); ylim([-1,1]*pi*doppler_to_speed)
            ylabel(['after ',num2str(i_target),' extraction(s)']);
        end
        figure(fig_0); hold on;
      
        plot((ind_col-1)*2*pi/N_symb/oversampling_symb*ffreq_to_range/ts,w_doppler ...
            *doppler_to_speed,'bo','LineWidth',1,'MarkerSize', ...
            5+5*power_bins(ind_row,ind_col)/peak_power)
    end
end
end