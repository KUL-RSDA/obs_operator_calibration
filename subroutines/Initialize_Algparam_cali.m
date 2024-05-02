%
% Set calibration algorithm parameters
%
% option: algorithm choice
%
% Last update:
% 12 Dec 2013: Gabrielle De Lannoy
% 23 Sep 2014: fed in select_par as input
% 20 Jun 2021: Michel Bechtold, extended w/ sceua
% 14 Mar 2023: Gabrielle De Lannoy, minor clean up, added subroutines
%==========================================================================

function [alg_param] = Initialize_Algparam_cali(...
    cali_option,select_par,path_subroutines)


%Specific for steepest gradient
%------------------------------
if (strcmp(cali_option,'steepest_gradient'))
    
    alg_param.dfac     = 0.0001;
    alg_param.psi      = 0.005;
    alg_param.N_iter   = 500;
    error([cali_option,' not implemented'])
    
%Specific for PSO
%------------------------------   
elseif (strcmp(cali_option,'pso'))
    
    alg_param.w0       = 0.9;  		% initial inertia weight
    alg_param.wt       = 0.7;  		% final inertia weight
    alg_param.c1       = 0.7;  		% cognitive parameter (influence best eigen position)
    alg_param.c2       = 1.3;  		% social parameter (influence best global position)
    alg_param.pop_size = 20;	   	% population size
    alg_param.delta    = 0.6; 		% fraction of par_size to determine vmax
    
    alg_param.N_iter   = 4;      	% minimum number of iterations (minimally 4!!)
    alg_param.N_iter_max = 20; 		% maximum number of iterations
    
    alg_param.OF_tol   = 1E-5; 		% after minim number of iterations, check if OF changes by more than OF_tol
    
    alg_param.repeat   = 4; 		% number of repetitions (~ensembles)
    
%SPecific for SCE-UA    
%------------------------------
elseif (strcmp(cali_option,'sceua'))
    
    % hard coded 5 parameters to adjust sceua algorithm parameters
    nopt=5;
    alg_param.ngs=1+2*nopt;
    alg_param.maxn=10000;
    alg_param.kstop=10;
    alg_param.pcento=0.1;
    alg_param.peps=0.001;
    alg_param.iniflg=0;
    rng('default'); % for same seed
    alg_param.seed = rng('shuffle'); % for new random numbers in each run
    
%Specific for MCMC
%------------------------------
elseif (strcmp(cali_option,'MCMC'))
    
    % Recommended parameter settings (runDREAM.m)
    alg_param.nCR = 3;                        % Crossover values used to generate proposals (geometric series)
    alg_param.DEpairs = 3;                    % Number of DEpairs
    alg_param.steps = 100;                    % Number of steps before output writing of iteration
    alg_param.eps = 0.05;                     % Random error for ergodicity
    alg_param.outlierTest = 'IQR_test';       % What kind of test to detect outlier chains?
    alg_param.pJumpRate_one = 0.2;            % Probability of selecting a jumprate of 1 --> jump between modes
    alg_param.pCR = 'Yes';                    % Adaptive tuning of crossover values (Yes or No)
    alg_param.Restart = 'No';                 % Restart run (Yes or No)?
    alg_param.modout = 'No';                  % Return model (function) simulations of samples (Yes or No)?
 					      % No! not set up yet to select multiple outputfields 
					      % as obtained from H-/V-/6angle SMOS
    
    % Problem specific parameter settings
    alg_param.n = length(select_par);         % Dimension of the problem (N_par)
    alg_param.seq = 10;                       % Number of Markov Chains / sequences
    alg_param.ndraw = 4000;                   % Maximum number of function evaluations
    alg_param.T = 1;                          % Each Tth sample is collected in the chains
    
    alg_param.prior = 'PRIOR';                % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    alg_param.prior_type = 'GAUSSIAN';        % Shape of prior parameter pdf: GAUSSIAN: "normrnd" or UNIFORM: "unifrnd")
    
    alg_param.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    alg_param.lik = 12;                       % 12: used for WCM calibration
                                              % 11: the model returns the multi-angular/pol TB OF
    
    % Provide information to do covariance sampling ("COV")
    % Currently pure placeholder and filled in later in dream_RTM.m
    alg_param.mu = zeros(1,alg_param.n);        % Provide mean of initial sample
    alg_param.cov = eye(alg_param.n);           % Initial covariance (identity matrix)
    
    
%Specific for DREAM_ZS
%------------------------------
elseif (strcmp(cali_option,'DREAMZS'))
    
    alg_param.k = 10;                         % Thinning parameter for appending X to Z
    alg_param.parallelUpdate = 0.9;           % Fraction of parallel direction updates
    alg_param.steps = inline('MCMCPar.ndraw/(20 * MCMCPar.seq)'); % Number of steps before calculating convergence diagnostics
    alg_param.m0 = inline('10 * MCMCPar.n');  % Initial size of matrix Z
    
    alg_param.nCR = 3;                        % Crossover values used to generate proposals (geometric series)
    alg_param.DEpairs = 3;                    % Number of DEpairs
    %MCMCPar.steps = 100;                     % Number of steps before output writing of iteration
    alg_param.eps = 0.05;                     % Random error for ergodicity
    alg_param.outlierTest = 'IQR_test';       % What kind of test to detect outlier chains?
    alg_param.pJumpRate_one = 0.2;            % Probability of selecting a jumprate of 1 --> jump between modes
    alg_param.pCR = 'Yes';                    % Adaptive tuning of crossover values (Yes or No)
    alg_param.Restart = 'No';                 % Restart run (Yes or No)?
    alg_param.modout = 'No';                  % Return model (function) simulations of samples (Yes or No)?
				              % No! not set up yet to select multiple outputfields 
                                              % as obtained from H-/V-/6angle SMOS
    alg_param.save = 'No';                    % Save output during the run (Yes or No)
    alg_param.ABC = 'No';                     % Approximate Bayesian Computation or Not?
    
    % Problem specific parameter settings
    alg_param.n = length(select_par);         % Dimension of the problem (N_par)
    alg_param.seq = 3;                        % Number of Markov Chains / sequences
    alg_param.ndraw = 12000;                  % Maximum number of function evaluations
    alg_param.T = 1;                          % Each Tth sample is collected in the chains
    
    %MCMCPar.prior = 'COV'; was 'COV' for initial RSE submission of De Lannoy et al. 2014,
    %but erroneous and corrected
    alg_param.prior = 'PRIOR';                % Latin Hypercube sampling (options, "LHS", "COV" and "PRIOR")
    alg_param.prior_type = 'GAUSSIAN';        % Shape of prior parameter pdf: GAUSSIAN: "normrnd" or UNIFORM: "unifrnd")
    
    alg_param.BoundHandling = 'Reflect';      % Boundary handling (options, "Reflect", "Bound", "Fold", and "None");
    alg_param.lik = 11;                       % The model returns the multi-angular/pol TB OF
    					      % (could be Tb_clim or Tb_MSE)   
    
    % Provide information to do covariance sampling ("COV")
    % Currently pure placeholder and filled in later in dream_RTM.m
    alg_param.mu    = zeros(1,alg_param.n);   % Provide mean of initial sample
    alg_param.cov   = eye(alg_param.n);       % Initial covariance (identity matrix)
    alg_param.N_ens = 20;                     % To save a subsample of ensemble output for uncertainty estimation
                                              % (param_ens)
    
    % Calculate MCMCPar.steps
    alg_param.steps = floor(alg_param.steps(alg_param));
    
    % Calculate MCMCPar.m0
    alg_param.m0 = alg_param.m0(alg_param);
    
elseif (strcmp(cali_option,'random_search') || ...
        strcmp(cali_option,'SCE_gent') || ...
        strcmp(cali_option,'nocali'))
    
    alg_param = NaN;
    
else
    
    error('this calibration option does not exist')
    
end

if (strcmp(cali_option,'sceua'))
    addpath([path_subroutines,'/SCEUA']); 
elseif (strcmp(cali_option,'pso'))
    addpath([path_subroutines,'/PSO']); 
elseif (strcmp(cali_option,'SCE_gent'))
    addpath([path_subroutines,'/SCE_gent']); 
elseif (strcmp(cali_option,'DREAMZS'))
    addpath([path_subroutines,'/DREAMZS']); 
    addpath([path_subroutines,'/DREAMZS/MATLAB-Code-DREAM_ZS-Sequential-V1.4']); 
end

end
%=============================EOF===============================


