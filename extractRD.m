%% Two Dinensioal NOMP functions

function [omegaRangeList,omegaDopplerList, gainList, residueList] = extractRD(y,y_RD,...
			      	   tau, numRefine, window, Nbins, Rmin, Rmax, doppler_to_speed, ...
                       range_axis, speed_axis, overSamplingRate_R, overSamplingRate_D)
% SUMMARY:
% 
%   given measurements: y = S * (mixture of sinusoids) + white noise
%          returns two **ordered** lists: a list of 
%          estimates of frequencies of the sinusoids in the mixture 
%          and list of corresponding gains 
% INPUT:
%    y - measurements (B x L x M) (N_beacons x N_chirp x N_symbols)
%    S - measurement matrix (N.A for Identity matix)
%    tau - minimum payoff per sinusoid
%        - should be slightly more than the noise level sigma2
%        minimize: norm(y - sum of sinusoids)^2 + tau * ell0_norm(gain)
%        ** LARGE VALUES OF TAU PROMOTE SPARSITY **% 
%    overSamplingRate (optional) -  # grid points = overSamplingRate * N
%              
%    numRefine (optional) : number of refinement steps after each frequency
%                           
% OUTPUT:
%   omegaList    - frequencies
%   gainList     - gains of estimated frequencies
%   residueList  - trajectory of the energy in the residual


if ~exist('numRefine','var'), numRefine = 6;
elseif isempty(numRefine),    numRefine = 6; end

if (~exist('opt', 'var'))
        opt = true;
end

M = size(y,3); % N_symbols
L = size(y,2); % N_chirps
N_beacons = size(y,1);


% Algorithm preprocessing
    sampledManifold = preProcessMeasMat(M,L, overSamplingRate_R, overSamplingRate_D ); 

if sampledManifold.is_eye
	S = []; % erase identity matrix
end

omegaRangeList = [];
omegaDopplerList = [];
gainList  = [];

y_r = y;
y_rd_residue = y_RD;

power_residue = reshape(sum(abs(y_rd_residue).^2,1),L*overSamplingRate_D, M*overSamplingRate_R);
residueList = [ y_r(:)' * y_r(:) ];
res_infoList = [];
fig_2 = figure();
title("Extracted Power in NOMP");
timer = 100;
subplot_i = 1;

while timer>0
    % keep detecting new sinusoids until power in residue 
    % becomes small; *** how small *** determined by *** tau ***

    

    % detect gain and frequency of an additional sinusoid
    [omega_new_range, omega_new_doppler, y_r, y_rd_residue, power_residue,...
        h_l, res_inf_normSq_rot] =...
        detectNewRD(y_r, y_rd_residue, sampledManifold, window, power_residue);
    % detecttNew removes the contribution of the newly detected
    % from the old residue  y_r(input) and reports the new residual measurement y_r (output)
    res_infoList = [res_infoList ; res_inf_normSq_rot];
    
    % residue
   subplot_i = plotResidue(power_residue, timer,  range_axis, speed_axis,...
                fig_2, Rmin, Rmax, doppler_to_speed, subplot_i);

    % newly detected sinusoid is coarse - so we refine it to 
    % imitate detection on the continuum
    for i = 1:numRefine
        [omega_new_range, omega_new_doppler, h_l, y_r] = refineOne(y_r, omega_new_range,...
            omega_new_doppler, h_l, sampledManifold.ant_idx_range,...
            sampledManifold.ant_idx_doppler, true,window);
    end
    % refineOne checks whether the refinement step decreases
    % the l-2 norm of the residue y_r

    % Add newly detected sinusoid to the ordered lists
    omegaRangeList = [omegaRangeList; omega_new_range];
    omegaDopplerList = [omegaDopplerList; omega_new_doppler];
    gainList  = [gainList; power_residue];

    % refine all frequencies detected so far
    % can be interpreted as a search for better frequency supports
    [omegaRangeList,omegaDopplerList, gainList, y_r] = refineAll(y_r, omegaRangeList,...
        omegaDopplerList, gainList, sampledManifold, numRefine,window);
    % refineAll only uses refineOne to tweak parameters and the energy 
    % in the residual measurements y_r can only decrease as a result

    % Solve least squares for the dictionary set [Ax(omega)] omega in 
    % omegaList
    % This is taken care by the refineAll function
%     [omegaList, gainList, y_r] = solveLeastSquares(y , omegaList, ...
%         S, sampledManifold.ant_idx);    
%     

    residue_new = y_r(:)'*y_r(:);
    residueList = [residueList; residue_new];
    timer = timer-1;

    % stopping criterion:
    if res_inf_normSq_rot < tau
        break;
    end
end

end


function [omega_range, omega_doppler, y_r, y_rd_residue, power_residue, h_l, res_inf_normSq_rot] = detectNewRD(y_r, y_rd_residue, ...
					 sampledManifold,window, power_residue)
% SUMMARY:
% 	detects a new sinusoid on the coarse grid
% INPUT:
% 	y_r - measurements (N_beacons x N_chirps, N_symbols)
% 	sampledManifold - precomputed IFFTs (if measurements are compressive)
%
% OUTPUT:
% 	omega - frequency on [0,2*pi) which best explains the measurements y_r
% 	h_l  - corresponding complex gain
% 	r_l   - after removing the detected sinusoid from the
%         	measurements y, r_l is the ***residual measurement***
%   res_inf_normSq_rot - max energy among DFT directions - needed
%          for stopping criterion
    R1 = length(sampledManifold.coarseOmega_range);
    R2 = length(sampledManifold.coarseOmega_doppler);
    N = sampledManifold.length; % [M,L]
    OSR = [round(R2/N(2)),round(R1/N(1))] ; 
    N_chirp = N(2);
    N_beacons = size(y_r,1);
    N_symb = N(1);
%     for i_beacon = 1:N_beacons
%        % TODO : verify   vvvvvvvvvvvvv
%         y_rd(i_beacon,:,:) = reshape(fft2(reshape(y(i_beacon,:,:),N_chirp,N_symb), ...
%                          R2,R1),1,R2,R1);
%     end
    % aggregate power in each range-doppler bin from all beacons
    
    peak_power = max(max(power_residue));
   
   [max_v,ind_row_vect] = max(power_residue);
   [~,ind_col] = max(max_v);
   ind_row = ind_row_vect(ind_col);
   % Detect the range and doppler values
   omega_range = sampledManifold.coarseOmega_range(ind_col);
   omega_doppler = angle(exp(-1j*sampledManifold.coarseOmega_doppler(ind_row)));
 

   % compute the response corresponding to the
   % current estimate of the sinusoid
        
   x = reshape(exp(-1i*omega_doppler*(0:N_chirp-1)')* ...
        exp(1i*omega_range*(0:N_symb-1)),1,N_chirp,N_symb);
    
   if (window)
        for chirp_i = 1:N_chirp
            x(1,chirp_i,:)  = reshape(x(1,chirp_i,:),N_symb,1).* hann(length(x)) ;   
        end
        %%% Window the signal in doppler domain
        for syms_i = 1:N_symb
            x(1,:,syms_i)= reshape(x(1,:,syms_i),N_chirp,1).*hann(size(x,2));
        end
   end

    for i_beacon = 1:N_beacons % extract the identified target's response 
        % from signal residue in fast-slow time and range-doppler domain
        % residual measurements after subtracting
         % out the newly detected sinusoid
        y_r(i_beacon,:,:) = y_r(i_beacon,:,:) - ...
            y_rd_residue(i_beacon,ind_row,ind_col)/N_symb/N_chirp*x;
         y_rd_residue(i_beacon,:,:) = reshape(fft2(reshape(y_r(i_beacon,:,:), ...
            N_chirp,N_symb),R2,R1),1,R2,R1);
    end

    power_residue = reshape(sum(abs(y_rd_residue).^2,1),R2, R1);
   
   h_l = sum(y_rd_residue(:,ind_row, ind_col),1)/N_symb/N_chirp;% ... TODO Verify
   h_l = sum(squeeze(x(:))' * y_r(:,:)' / norm(squeeze(x))^2);
   % For stopping criterion
   % we check only DFT frequencies - gives us a handle
   % the false alarm rate (a measure of how often we would
   % over estimate the size of the support set)
  res_inf_normSq_rot = sqrt(max(max(power_residue(1:OSR(1):end, 1:OSR(2):end)))); %TODO: verify this


end

% --------------------------------------------------------------------

function [omega_r,omega_d, h_l, y_r] = refineOne(y_r, omega_r, omega_d, h_l,...
			 ant_idx_r, ant_idx_d, isOrth,window)
% SUMMARY:
%   Refines parameters (gain and frequency) of a single sinusoid
%   and updates the residual measurement vector y_r to reflect
%   the refinement
% INPUT:
% 	y_r - residual measurement (all detected sinusoids removed)
%	omega - current estimate of frequency of sinusoid we want to
% 			refine
%	gain - current estimate of gain of sinusoid we want to refine
% 	ant_idx - translating indexes to phases in definition of sinusoid
%	isOrth - binary flag - is y_r orthogonal to x(omega) 
%	       - default - false
% OUTPUT:
%       refined versions of omega, gain and y_r
%       (see INPUT for definitions)

    if ~exist('isOrth', 'var'), isOrth = false; end
      
    
    M = length(ant_idx_r);
    L = length(ant_idx_d);
    N_beacons = size(y_r,1);

    x_theta = reshape(exp(-1i*omega_d*(0:L-1)')* ...
        exp(1i*omega_r*(0:M-1)),1,L,M);%  TODO:Verify
    
   if (window)
        for chirp_i = 1:L
            x_theta(1,chirp_i,:)  = reshape(x_theta(1,chirp_i,:),M,1).* hann(length(x_theta)) ;   
        end
        %%% Window the signal in doppler domain
        for syms_i = 1:M
            x_theta(1,:,syms_i)= reshape(x_theta(1,:,syms_i),L,1).*hann(size(x_theta,2));
        end
   end

    % Differentiate wrt range and doppler
    % dx_theta = [d(x_theta)/d(omega_r); d(x_theta)/d(omega_d)]
    dx_theta_domega_r = (1j * squeeze(x_theta) .* ant_idx_r'); % (L x M)
    dx_theta_domega_d = (1j * ant_idx_d .* squeeze(x_theta)); % (Lx M)
    
    % Differentiate twice 
    % d2x_theta = [d^2(x_theta)/d^2(omega_r);d^2(x_theta)/d(omega_r)/d(omega_d);
    %              d^2(x_theta)/d(omega_r)/d(omega_d); d^2(x_theta)/d^2(omega_d)]
    d2x_theta_d2omega_r = (-1j * squeeze(x_theta) .* (ant_idx_r.^2)');
    %d2x_theta_domega_r_domega_d = -1j * ant_idx_d .* x_theta .*ant_idx_r';
    d2x_theta_d2omega_d = -1j * (ant_idx_d.^2) .* squeeze(x_theta);

%     d2x_theta =  cat(3,d2x_theta_d2omega_r(:)', d2x_theta_domega_r_domega_d(:)',...
%                   d2x_theta_domega_r_domega_d(:)' , d2x_theta_d2omega_d(:)'); % (LM x 1 x 4)
    
    % add the current estimate of the sinusoid to residue
    y = y_r + h_l.*x_theta; 
    
    % UPDATE GAIN
    % recompute gain and residue to ensure that 
    % y_r is orthogonal to x_theta - this property is lost when we refine other sinusoids
    if ~isOrth
        %gain = (x_theta'*y);
        y_r = y - h_l.*x_theta;
    end
    
    for i_beacon = 1:N_beacons
        der1_r = -2*real(h_l .* y_r(i_beacon,:)*dx_theta_domega_r(:)); 
        der1_d = -2*real(h_l .* y_r(i_beacon,:)*dx_theta_domega_d(:)); 
    
        der2_r = -2*real(h_l * y_r(i_beacon,:)*d2x_theta_d2omega_r(:)) +...
        2*(abs(h_l)^2)*(dx_theta_domega_r(:)'*dx_theta_domega_r(:));
    
        der2_d = -2*real(h_l * y_r(i_beacon,:)*d2x_theta_d2omega_d(:)) +...
        2*(abs(h_l)^2)*(dx_theta_domega_d(:)'*dx_theta_domega_d(:));
    end
    
    % UPDATE Range and doppler estiamtes 
   
    if der2_r > 0
        omega_r_next = omega_r - der1_r/der2_r;
    else
        omega_r_next = omega_r - sign(der1_r)*(1/4)*(2*pi/M)*rand(1);
        disp('error :second derivative of range is negative');
    end

    if der2_d > 0
        omega_d_next = omega_d - der1_d/der2_d;
    else
        omega_d_next = omega_d - sign(der1_d)*(1/4)*(2*pi/L)*rand(1);
        disp('error :second derivative of doppler is negative');
    end

    % COMPUTE x_theta for omega_next so that we can compute 
    % gains_next and y_r_next
    x_theta = reshape(exp(-1i*omega_d_next*(0:L-1)')* ...
        exp(1i*omega_r_next*(0:M-1)),1,L,M);
    
    if (window)
        for chirp_i = 1:L
            x_theta(1,chirp_i,:)  = reshape(x_theta(1,chirp_i,:),M,1).* hann(length(x_theta)) ;   
        end
        %%% Window the signal in doppler domain
        for syms_i = 1:M
            x_theta(1,:,syms_i)= reshape(x_theta(1,:,syms_i),L,1).*hann(size(x_theta,2));
        end
    end
    
    for i_beacon = 1:N_beacons
        % UPDATE GAIN
        gain_next = (x_theta(:)'*y(i_beacon,:)');
        
        % UPDATE RESIDUE
        y_r_next(i_beacon,:,:) = y(i_beacon) - gain_next*x_theta;
        
        % check for decrease in residue -  needed as a result of 
        % non-convexity of residue (even when the cost surface 
        % is locally convex); This is the same as checking whether 
        % |<y, x_theta>|^2/<x_theta,x_theta> improves as we ensured 
        % that y - gain*x_theta is perp to x_theta by recomputing gain
        if (y_r_next(i_beacon,:)*y_r_next(i_beacon,:)') <= (y_r(i_beacon,:)*y_r(i_beacon,:)')
            % commit only if the residue decreases
            omega_r = omega_r_next;
            omega_d = omega_d_next;
            h_l = gain_next;
            y_r(i_beacon) = y_r_next(i_beacon);
        end
    end
end

% --------------------------------------------------------------------

function [omegaRangeList,omegaDopplerList, gainList, y_r] = refineAll(y_r, omegaRangeList,...
        omegaDopplerList, gainList, sampledManifold, numRefine, window)
% SUMMARY:
%   uses refineOne algorithm to refine frequencies and gains of
%   of all sinusoids
% INPUT:
%    y_r - residual measurement after all detected sinusoids have been
%          removed
%    omegaRangeList - list of range frequencies of known(detected) sinusoids
%    omegaDopplerList - list of doppler frequencies of known(detected) sinusoids
%    gainList  - list of gains of known(detected) sinusoids
%  
%    ant_idx - translating indexes to phases in definition of sinusoid
%    numRefine - number of times each sinusoid in the list is refined
%
% OUTPUT:
%       refined versions of inputs omegaList, gainList, y_r

K = length(omegaRangeList); % number of sinusoids

% refinement repeated "numRefine" number of times per sinusoid

for i = 1:numRefine
    % chose an ordering for refinement
    
    % % *random* ordering
    % order = randperm(K);
    
    % *sequential* ordering
    order = 1:K;
    
    for j = 1:K
        l = order(j);
        
        % parameters of the l-th sinusoid
        omega_r = omegaRangeList(l);
        omega_d = omegaDopplerList(l);
        gain = gainList(l);
        
        % refine our current estimates of (gain, omega) of the
        % l-th sinusoid
        [omega_r,omega_d, gain, y_r] = refineOne(y_r, omega_r, omega_d, gain,...
			 sampledManifold.ant_idx_range, sampledManifold.ant_idx_doppler, false, window);
      
        omegaRangeList(l) = omega_r;
        omegaDopplerList(l) = omega_d;
        gainList(l) = gain;
        % refineOne ensures that (gain, omega) pair are altered iff
        % the energy in the residual measurements y_r does not 
    	% increase
    end
    
end

end

% --------------------------------------------------------------------

% --------------------------------------------------------------------

function sampledManifold = preProcessMeasMat(M,L, overSamplingRate_R, overSamplingRate_D)
% ***WIP : Use for only non-compressive measurements****
% SUMMARY:
%   compute overSamplingRate*(M,L) IFFTs once 
% INPUT:
%       M - length of sinusoids in range axis
%       L - length of sinusoids in doppler axis
%    overSamplingRate (optional) - how fine (in terms of multiples
%       of the FFT grid) do we want the coarse grid to be?
%       number of grid points = overSamplingRate * M * L
% OUTPUT:
%       data structure with sinusoid responses
%       precomputed using IFFTs 
    
    sampledManifold.length = [M,L];
    R1 = round(overSamplingRate_R*M);
    R2 = round(overSamplingRate_D*L);
    
    sampledManifold.oversampled_range = R1;
    sampledManifold.oversampled_doppler = R2;
    sampledManifold.coarseOmega_range = 2*pi*(0:(R1-1))/R1 ;
    sampledManifold.coarseOmega_doppler = 2*pi*(0:(R2-1))/R2 ;  % omegaCoarse
    
    % definition of a sinusoid:
    %              exp(-1j*omega((N-1)/2):((N-1)/2))/sqrt(N)
    
    ant_idx_range = 0:(M-1);
    ant_idx_doppler = 0:(L-1);
    
    % IFFT definition of a sinusoid(omega) takes the following form:
    % 	sinusoid    = @(omega) exp(1j*(0:(M-1)).'*omega);
    % To reiterate, we assume that a sinusoid is given by
    %	sinusoid    = @(omega) exp(1j*ant_idx.'*omega)/sqrt(N);
    % So we store this information in sampledManifold container
    sampledManifold.ant_idx_range = ant_idx_range(:);
    sampledManifold.ant_idx_doppler = ant_idx_doppler(:);
    
    sampledManifold.is_eye = true;

end

function subplot_i = plotResidue(power_residue, iteration , range_axis,speed_axis ,fig_2,...
    Rmin, Rmax, doppler_to_speed, subplot_i)
    
    if mod(iteration, 4) == 0
        figure(fig_2); hold on;
        subplot(2,2,subplot_i);
        h = surf(range_axis,speed_axis,10*log10(power_residue));
        set(h,'edgecolor','none');view(3);
        % xlim([0,N_symb*oversampling_symb]); ylim([0,N_chirp*oversampling_chirp]);
        xlim([Rmin,Rmax]); ylim([-1,1]*pi*doppler_to_speed)
        ylabel(['after ',num2str(iteration),' extraction(s)']);
        subplot_i = subplot_i +1;
    end
end
