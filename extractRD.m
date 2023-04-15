%% Two Dinensioal NOMP functions

function [omegaRangeList,omegaDopplerList, gainList, residueList] = extractRD(y,...
			      	   tau, overSamplingRate, numRefine)
% SUMMARY:
% 
%   given measurements: y = S * (mixture of sinusoids) + white noise
%          returns two **ordered** lists: a list of 
%          estimates of frequencies of the sinusoids in the mixture 
%          and list of corresponding gains 
% INPUT:
%    y - measurements (L x M) (N_chirp x N_symbols)
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

if ~exist('overSamplingRate','var'), overSamplingRate = 4;
elseif isempty(overSamplingRate), overSamplingRate = 4; end

if ~exist('numRefine','var'), numRefine = 6;
elseif isempty(numRefine),    numRefine = 6; end

M = size(y,2); % N_symbols
L = size(y,1); % N_chirps
overSamplingRate_R = 8;
overSamplingRate_D = 32;
% Algorithm preprocessing
    sampledManifold = preProcessMeasMat(M,L, overSamplingRate_R, overSamplingRate_D ); 

if sampledManifold.is_eye
	S = []; % erase identity matrix
end

omegaRangeList = [];
omegaDopplerList = [];
gainList  = [];
y_r = y;
residueList = [ y_r(:)' * y_r(:) ];
res_infoList = [];
timer = 10;
while timer>0
    % keep detecting new sinusoids until power in residue 
    % becomes small; *** how small *** determined by *** tau ***
    
    % detect gain and frequency of an additional sinusoid
    [omega_new_range, omega_new_doppler, h_l, y_r, res_inf_normSq_rot] =...
        detectNewRD(y_r, sampledManifold);
    % detecttNew removes the contribution of the newly detected
    % from the old residue  y_r(input) and reports the new residual measurement y_r (output)
    res_infoList = [res_infoList ; res_inf_normSq_rot];
    % stopping criterion:
    if res_inf_normSq_rot < tau
        break;
    end

    % newly detected sinusoid is coarse - so we refine it to 
    % imitate detection on the continuum
    for i = 1:numRefine
        [omega_new_range, omega_new_doppler, h_l, y_r] = refineOne(y_r, omega_new_range,...
            omega_new_doppler, h_l, sampledManifold.ant_idx_range,...
            sampledManifold.ant_idx_doppler, true);
    end
    % refineOne checks whether the refinement step decreases
    % the l-2 norm of the residue y_r

    % Add newly detected sinusoid to the ordered lists
    omegaRangeList = [omegaRangeList; omega_new_range];
    omegaDopplerList = [omegaDopplerList; omega_new_doppler];
    gainList  = [gainList; h_l];

    % refine all frequencies detected so far
    % can be interpreted as a search for better frequency supports
    [omegaRangeList,omegaDopplerList, gainList, y_r] = refineAll(y_r, omegaRangeList,...
        omegaDopplerList, gainList, sampledManifold, numRefine);
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
end

end


function [omega_range, omega_doppler, h_l, r_l, res_inf_normSq_rot] = detectNewRD(y,...
					 sampledManifold)
% SUMMARY:
% 	detects a new sinusoid on the coarse grid
% 	refines the newly detected frequency *numRefine* times
% 
% INPUT:
% 	y - measurements (N_chirps, N_symbols)
% 	sampledManifold - precomputed IFFTs (if measurements are compressive)
%
% OUTPUT:
% 	omega - frequency on [0,2*pi) which best explains the measurements y
% 	h_l  - corresponding complex gain
% 	r_l   - after removing the detected sinusoid from the
%         	measurements y, r_l is the ***residual measurement***
%   res_inf_normSq_rot - max energy among DFT directions - needed
%          for stopping criterion
    R1 = length(sampledManifold.coarseOmega_range);
    R2 = length(sampledManifold.coarseOmega_doppler);
    N = sampledManifold.length; % [M,L]
    OSR = [round(R2/N(2)),round(R1/N(1))] ; 
  
   gains  = fft2(y, R2, R1)/sqrt(N(1))/sqrt(N(2)); % Non compressive measurements do not need
                % <y x(w)> / ||x(w)|| 
   prob = (abs(gains).^2);
  
   [max_v,ind_row_vect] = max(prob);
   [~,ind_col] = max(max_v);
   ind_row = ind_row_vect(ind_col);

   omega_range = sampledManifold.coarseOmega_range(ind_col);
   omega_doppler = angle(exp(-1j*sampledManifold.coarseOmega_doppler(ind_row))); % TODO : Verify
 
   gain = gains(ind_col,ind_row);

   % compute the response corresponding to the
   % current estimate of the sinusoid
        
   x = exp(-1j*sampledManifold.ant_idx_doppler * omega_doppler)*...
       exp(1j*sampledManifold.ant_idx_range * omega_range)'/...
       sqrt(sampledManifold.length(1))/sqrt(sampledManifold.length(2)); %(N-chirp x N_symbol) TODO:Verify
    
   h_l = x(:)' * y(:) / norm(x)^2; 
   % residual measurements after subtracting
   % out the newly detected sinusoid
   r_l  = y - h_l .* x;

   % For stopping criterion
   % we check only DFT frequencies - gives us a handle
   % the false alarm rate (a measure of how often we would
   % over estimate the size of the support set)

   res_inf_normSq_rot = sqrt(max(max(prob(1:OSR(1):end, 1:OSR(2):end)))); %TODO: verify this
end

% --------------------------------------------------------------------

function [omega_r,omega_d, h_l, y_r] = refineOne(y_r, omega_r, omega_d, h_l,...
			 ant_idx_r, ant_idx_d, isOrth)
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
    x_theta  = exp(1j*ant_idx_d*omega_d)*exp(1j*ant_idx_r*omega_r)'/sqrt(M)/sqrt(L); % (L x M)
    % Differentiate wrt range and doppler
    % dx_theta = [d(x_theta)/d(omega_r); d(x_theta)/d(omega_d)]
    dx_theta_domega_r = (1j * x_theta .* ant_idx_r'); % (L x M)
    dx_theta_domega_d = (1j * ant_idx_d .* x_theta); % (Lx M)
    
    % Differentiate twice 
    % d2x_theta = [d^2(x_theta)/d^2(omega_r);d^2(x_theta)/d(omega_r)/d(omega_d);
    %              d^2(x_theta)/d(omega_r)/d(omega_d); d^2(x_theta)/d^2(omega_d)]
    d2x_theta_d2omega_r = (-1j * x_theta .* (ant_idx_r.^2)');
    %d2x_theta_domega_r_domega_d = -1j * ant_idx_d .* x_theta .*ant_idx_r';
    d2x_theta_d2omega_d = -1j * (ant_idx_d.^2) .* x_theta;

%     d2x_theta =  cat(3,d2x_theta_d2omega_r(:)', d2x_theta_domega_r_domega_d(:)',...
%                   d2x_theta_domega_r_domega_d(:)' , d2x_theta_d2omega_d(:)'); % (LM x 1 x 4)
    
    % add the current estimate of the sinusoid to residue
    y = y_r + h_l.*x_theta; 
    
    % UPDATE GAIN
    % recompute gain and residue to ensure that 
    % y_r is orthogonal to x_theta - this property is lost when we refine other sinusoids
    if ~isOrth
        gain = (x_theta'*y);
        y_r = y - h_l.*x_theta;
    end
    
    der1_r = -2*real(h_l .* y_r(:)'*dx_theta_domega_r(:)); 
    der1_d = -2*real(h_l .* y_r(:)'*dx_theta_domega_d(:)); 

    der2_r = -2*real(h_l * y_r(:)'*d2x_theta_d2omega_r(:)) +...
    2*(abs(h_l)^2)*(dx_theta_domega_r(:)'*dx_theta_domega_r(:));

    der2_d = -2*real(h_l * y_r(:)'*d2x_theta_d2omega_d(:)) +...
    2*(abs(h_l)^2)*(dx_theta_domega_d(:)'*dx_theta_domega_d(:));

    
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
    x_theta  = exp(1j*ant_idx_d*omega_d_next)*...
        exp(1j*ant_idx_r*omega_r_next)'/sqrt(M)/sqrt(L); % (L x M)

    % UPDATE GAIN
    gain_next = (x_theta(:)'*y(:));
    
    % UPDATE RESIDUE
    y_r_next = y - gain_next*x_theta;
    
    % check for decrease in residue -  needed as a result of 
    % non-convexity of residue (even when the cost surface 
    % is locally convex); This is the same as checking whether 
    % |<y, x_theta>|^2/<x_theta,x_theta> improves as we ensured 
    % that y - gain*x_theta is perp to x_theta by recomputing gain
    if (y_r_next(:)'*y_r_next(:)) <= (y_r(:)'*y_r(:))
        % commit only if the residue decreases
        omega_r = omega_r_next;
        omega_d = omega_d_next;
        h_l = gain_next;
        y_r = y_r_next;
    end

end

% --------------------------------------------------------------------

function [omegaRangeList,omegaDopplerList, gainList, y_r] = refineAll(y_r, omegaRangeList,...
        omegaDopplerList, gainList, sampledManifold, numRefine)
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
			 sampledManifold.ant_idx_range, sampledManifold.ant_idx_doppler, false);
      
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