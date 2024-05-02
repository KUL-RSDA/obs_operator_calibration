%% Particle Swarm Optimization (PSO)
%
% Gabrielle De Lannoy, 11 March 2011
% Gabrielle De Lannoy, 25 November 2019, 
%                      largely rewritten & generalized
% Sara Modanesi,       25 August 2020 modified for SSE
% Gabrielle De Lannoy, 13 March 2023, clean up
%
% OUTPUT
%---------------------------------------------------------
% param_t_0   -- output: structure w/ optimal parameter values
% OF_iter_t_0 -- output: value of OF_t for all iterations (last = optimal)
%                        (total OF)
% OF_iter_m_0 -- output: value of OF_m for all iterations (last = optimal)
%                        (sum of square error as part of OF)
% OF_iter_p_0 -- output: (parameter penalty term)   
% 
% INPUT
%--------------------------------------------------------
% pso         -- input: structure w/ PSO-specific parameters
% fwd_in      -- input: structure with 
%                       - time series of dynamic fields for a single pixel (e.g. sm, lai, backscatter [obs],...)
% obs         -- input: array or structure of obs (in shape expected by OF)
%                       - e.g. Sentinel-1 backscatter
% param_0     -- input: structure w/ initial value (m) and range (r) 
%                       for all parameters (used as prior constraint)  
% select_par  -- input: array with indices selecting parameters for
%                       calibration relative to the total array of param_0
%                       (control parameters)
% OF_param    -- input: structure of constants needed in OF computation
%                       (e.g. nopriorscaling, climscaling, ...
%                       W_m, W_s, W_p, ...
%                       name it to what you need in the OF of your choice)
%
%=========================================================
                
function [param_t_0, OF_iter_t_0, OF_iter_m_0, OF_iter_p_0 ] = ...
    pso_WCM_SSE(pso, ...
    fwd_in, ...
    obs, ...
    param_0, select_par, ...    
    OF_param)

%Algorithm initialization

OF_t_min    = 999999999999;

fn          = fieldnames(param_0.m);
par_size    = length(select_par);               % number of parameters

x_          = zeros([pso.pop_size, par_size]);  % instantaneous particle position
OF_t        = zeros(pso.pop_size,3);            % corresponding OF (composed by three terms: total, SSE, param penalty)
p_          = zeros([pso.pop_size,par_size]);   % best particle position till now
pOF_t       = zeros(pso.pop_size,1);            % corresponding OF (total ONLY!)

pOF_sse     = zeros(pso.pop_size,1);             % corresponding OF (SSE) 
pOF_p       = zeros(pso.pop_size,1);            % corresponding OF (p)  

v_          = zeros([pso.pop_size, par_size]);  % instantaneous particle speed
range       = zeros(par_size,1);

p_constraint.var     = param_0.m;  %Initialize prior variance
p_constraint.pen     = param_0.m;  %Constraining parameter value, i.e. "prior"

rand('state',sum(100*clock)); % initialize RN generator

if ~(strcmp(OF_param.OF_type,'SSE'))
    error(['not ready for OF-option ',OF_param.OF_type])    
end

%Loop through repetitions

for l=1:pso.repeat

    param_0_tmp = param_0.m;
    OF_iter     = NaN+zeros(pso.N_iter_max+1,3); %3 components of OF
    par_cov     = NaN+zeros(pso.N_iter_max+1,par_size,par_size);

    %%Initialization

    %Initial position and velocity

    for p=1:par_size

        range(p)  = max([param_0.r.(fn{select_par(p)})]) - min([param_0.r.(fn{select_par(p)})]);  
        
        %starting from best guess may squeeze the ensemble in one corner and
        %slow down the convergence => start from the middle of the possible
        %range now instead

        %rand generates elements in [0 1]
        x_(:,p)   = rand(pso.pop_size,1)*range(p) + min([param_0.r.(fn{select_par(p)})]);            
                    
        v_(:,p)   = randn(pso.pop_size,1)*range(p) / pso.N_iter; %division by N_iter to avoid initial high velocity

        %reset outliers
        x_((x_(:,p) > max([param_0.r.(fn{select_par(p)})])),p) = max([param_0.r.(fn{select_par(p)})]);  % boundary cond
        x_((x_(:,p) < min([param_0.r.(fn{select_par(p)})])),p) = min([param_0.r.(fn{select_par(p)})]);  % boundary cond

        if (OF_param.W_p ~= 0) 
            % Prepare to add a parameter penalty
            p_constraint.var.(fn{select_par(p)}) = range(p).^2./12;  %variance of uniform distribution
            p_constraint.pen.(fn{select_par(p)}) = param_0.m.(fn{select_par(p)});
        end

    end

    %Initial OF

    for i=1:pso.pop_size

      for p=1:par_size  
        param_0_tmp.(fn{select_par(p)}) = x_(i,p);
      end

      % OF
      %--------------------
      [OF_tmp, OF_sse_tmp, OF_p_tmp] = ...
          CalWCM_SSE(fwd_in, ...
            param_0_tmp, obs, p_constraint, ...        
            OF_param, select_par);

      OF_t(i,1) = OF_tmp;
      OF_t(i,2) = OF_sse_tmp;
      OF_t(i,3) = OF_p_tmp;
      
    end

    p_        = x_;             % best position till now = start position
    pOF_t     = OF_t(:,1);      % corresponding best OF = start RMSE, ...
    pOF_sse   = OF_t(:,2);
    pOF_p     = OF_t(:,3);

    [pgOF_t,g] = min(pOF_t);     % index g for global best OF

    OF_iter(1,1) = pOF_t(g);
    OF_iter(1,2) = pOF_sse(g);
    OF_iter(1,3) = pOF_p(g);
    
    %Covariance (uncertainty) matrix 
    par_cov(1,:,:) = cov(p_);
        
    %%Iteration: adjust velocity and position
    for iter=1:pso.N_iter_max

        %Adjust inertia weight for each iteration
        w = (pso.w0-pso.wt)*(pso.N_iter_max-iter)/pso.N_iter_max + pso.wt; 

        for i=1:pso.pop_size

            for p=1:par_size

                r1 = rand;
                r2 = rand;

                v_(i,p) = w*v_(i,p)+ ...
                    pso.c1*r1*(p_(i,p)-x_(i,p))+ ...
                    pso.c2*r2*(p_(g,p)-x_(i,p));

                if (abs(v_(i,p)) > pso.delta*range(p))
                    v_(i,p) = sign(v_(i,p))*pso.delta*range(p);
                end
                x_(i,p) = x_(i,p)+v_(i,p);

                if (x_(i,p) > max([param_0.r.(fn{select_par(p)})]))
                    x_(i,p) = max([param_0.r.(fn{select_par(p)})]);  % boundary cond
                    v_(i,p) = -v_(i,p);      % velocity reflection
                end
                if (x_(i,p) < min([param_0.r.(fn{select_par(p)})]))
                    x_(i,p) = min([param_0.r.(fn{select_par(p)})]);  % boundary cond
                    v_(i,p) = -v_(i,p);      % velocity reflection
                end

                param_0_tmp.(fn{select_par(p)}) = x_(i,p);

            end

            %New OF
            %--------------------
            [OF_tmp, OF_sse_tmp, OF_p_tmp] = ...
            CalWCM_SSE(fwd_in, ...
            param_0_tmp, obs, p_constraint, ...        
            OF_param, select_par);
    
            OF_t(i,1) = OF_tmp;
            OF_t(i,2) = OF_sse_tmp;
            OF_t(i,3) = OF_p_tmp;

            if (OF_t(i,1) < pOF_t(i))
                p_(i,:)     = x_(i,:);    % update best position
                pOF_t(i)    = OF_t(i,1);

                pOF_sse(i) = OF_t(i,2);
                pOF_p(i)    = OF_t(i,3);
            end


        end

        [pgOF_t,g]      = min(pOF_t);  % index g of global best OF_t   

        OF_iter(iter+1,1) = pOF_t(g);

        OF_iter(iter+1,2) = pOF_sse(g);
        OF_iter(iter+1,3) = pOF_p(g);

        for p=1:par_size
          param_0_tmp.(fn{select_par(p)}) = p_(g,p);
        end

        %Covariance (uncertainty) matrix; population covar
        par_cov(iter+1,:,:) = cov(p_);
        
        %Get out of iterations if at bottom of OF_t
        if (iter >= pso.N_iter)

            if( abs(OF_iter(iter+1,1)-OF_iter(iter+1-3,1)) < pso.OF_tol)
                break
            end

        end

    end
    
    %Keep best parameters, only if this iteration is better than any previous one
    if (OF_iter(iter+1,1) < OF_t_min)

        OF_t_min = OF_iter(iter+1,1);

        param_t_0   = param_0_tmp;
        OF_iter_t_0 = OF_iter(:,1);
        OF_iter_m_0 = OF_iter(:,2);
        OF_iter_p_0 = OF_iter(:,3);
        
    end

end

%==========================EOF=======================================================
