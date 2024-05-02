%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [param_opt,bestf] = sceua_WCM(param_0,select_par,maxn,kstop,...
    pcento,peps,ngs,iseed,iniflg,fwd_in,obs,s_0,OF_param,debug)
% This is the subroutine implementing the SCE algorithm,
% written by Q.Duan, 9/2004
%
%  Michel Bechtold 20 Jun 2021, adapted for use with WCM
%  Gabrielle De Lannoy 14 Mar 2023, 
%       updated calls to OF, 
%       cast bestx to param_opt,
%       silenced reports (diplay only shown in debug mode)
%
% Definition:
%  x0 = the initial parameter array at the start;
%     = the optimized parameter array at the end;
%  f0 = the objective function value corresponding to the initial parameters
%     = the objective function value corresponding to the optimized parameters
%  bl = the lower bound of the parameters
%  bu = the upper bound of the parameters
%  iseed = the random seed number (for repetetive testing purpose)
%  iniflg = flag for initial parameter array (=1, included it in initial
%           population; otherwise, not included)
%  ngs = number of complexes (sub-populations)
%  npg = number of members in a complex
%  nps = number of members in a simplex
%  nspl = number of evolution steps for each complex before shuffling
%  mings = minimum number of complexes required during the optimization process
%  maxn = maximum number of function evaluations allowed during optimization
%  kstop = maximum number of evolution loops before convergency
%  percento = the percentage change allowed in kstop loops before convergency
% LIST OF LOCAL VARIABLES
%    x(.,.) = coordinates of points in the population
%    xf(.) = function values of x(.,.)
%    xx(.) = coordinates of a single point in x
%    cx(.,.) = coordinates of points in a complex
%    cf(.) = function values of cx(.,.)
%    s(.,.) = coordinates of points in the current simplex
%    sf(.) = function values of s(.,.)
%    bestx(.) = best point at current shuffling loop
%    bestf = function value of bestx(.)
%    worstx(.) = worst point at current shuffling loop
%    worstf = function value of worstx(.)
%    xnstd(.) = standard deviation of parameters in the population
%    gnrng = normalized geometri%mean of parameter ranges
%    lcs(.) = indices locating position of s(.,.) in x(.,.)
%    bound(.) = bound on ith variable being optimized
%    ngs1 = number of complexes in current population
%    ngs2 = number of complexes in last population
%    iseed1 = current random seed
%    criter(.) = vector containing the best criterion values of the last
%                10 shuffling loops
global BESTX BESTF ICALL PX PF

fn          = fieldnames(param_0.m);
par_size    = length(select_par);               % number of parameters
for p=1:par_size
    x0(p) = param_0.m.(fn{select_par(p)});
    bl(p) = param_0.r.(fn{select_par(p)})(1);
    bu(p) = param_0.r.(fn{select_par(p)})(2);
    % Prepare to add a parameter penalty
    p_constraint.var.(fn{select_par(p)}) = (bu(p)-bl(p)).^2./12;  %variance of uniform distribution
    p_constraint.pen.(fn{select_par(p)}) = param_0.m.(fn{select_par(p)});  
end

% Initialize SCE parameters:
nopt=length(x0);
npg=2*nopt+1;
nps=nopt+1;
nspl=npg;
mings=ngs;
npt=npg*ngs;
bound = bu-bl;

%Algorithm initialization

par_size    = length(x0);               % number of parameters

%%%
% Create an initial population to fill array x(npt,nopt):
rand('seed',double(iseed.State(1)));
x=zeros(npt,nopt);
for i=1:npt;
    x(i,:)=bl+rand(1,nopt).*bound;
end;
if iniflg==1; x(1,:)=x0; end;
nloop=0;
icall=0;

for i=1:npt;
    for p=1:par_size
        param_0.m.(fn{select_par(p)})  = x0(p);
    end
    if (strcmp(OF_param.OF_type,'SSE'))
        [OF_tmp, OF_sse_tmp, OF_p_tmp] = ...
            CalWCM_SSE(fwd_in, ...
            param_0.m, obs, p_constraint, ...
            OF_param,select_par);
    else
        error(['not ready for OF-option ',OF_param.OF_type])
    end
    xf(i) = OF_tmp;
    icall = icall + 1;
end;
f0=xf(1);
% Sort the population in order of increasing function values;
[xf,idx]=sort(xf);
x=x(idx,:);
% Record the best and worst points;
bestx=x(1,:); bestf=xf(1);
worstx=x(npt,:); worstf=xf(npt);
BESTF=bestf; BESTX=bestx;ICALL=icall;
% Compute the standard deviation for each parameter
xnstd=std(x);
% Computes the normalized geometric range of the parameters
gnrng=exp(mean(log((max(x)-min(x))./bound)));
if debug
disp('The Initial Loop: 0');
disp(['BESTF  : ' num2str(bestf)]);
disp(['BESTX  : ' num2str(bestx)]);
disp(['WORSTF : ' num2str(worstf)]);
disp(['WORSTX : ' num2str(worstx)]);
disp(' ');
end
% Check for convergency;
if debug
if icall >= maxn;
    disp('*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE LIMIT');
    disp('ON THE MAXIMUM NUMBER OF TRIALS ');
    disp(maxn);
    disp('HAS BEEN EXCEEDED.  SEARCH WAS STOPPED AT TRIAL NUMBER:');
    disp(icall);
    disp('OF THE INITIAL LOOP!');
end;
if gnrng < peps;
    disp('THE POPULATION HAS CONVERGED TO A PRESPECIFIED SMALL PARAMETER SPACE');
end;
end
% Begin evolution loops:
nloop = 0;
criter=[];
criter_change=1e+5;
while icall<maxn & gnrng>peps & criter_change>pcento;
    nloop=nloop+1;
    
    % Loop on complexes (sub-populations);
    for igs = 1: ngs;
        
        % Partition the population into complexes (sub-populations);
        k1=1:npg;
        k2=(k1-1)*ngs+igs;
        cx(k1,:) = x(k2,:);
        cf(k1) = xf(k2);
        
        % Evolve sub-population igs for nspl steps:
        for loop=1:nspl;
            
            % Select simplex by sampling the complex according to a linear
            % probability distribution
            lcs(1) = 1;
            for k3=2:nps;
                for iter=1:1000;
                    lpos = 1 + floor(npg+0.5-sqrt((npg+0.5)^2 - npg*(npg+1)*rand));
                    idx=find(lcs(1:k3-1)==lpos); if isempty(idx); break; end;
                end;
                lcs(k3) = lpos;
            end;
            lcs=sort(lcs);
            % Construct the simplex:
            s = zeros(nps,nopt);
            s=cx(lcs,:); sf = cf(lcs);
            
            [snew,fnew,icall]=cceua_WCM(s,sf,bl,bu,icall,maxn,fwd_in,obs,s_0,p_constraint,OF_param,param_0,select_par);
            % Replace the worst point in Simplex with the new point:
            s(nps,:) = snew; sf(nps) = fnew;
            
            % Replace the simplex into the complex;
            cx(lcs,:) = s;
            cf(lcs) = sf;
            
            % Sort the complex;
            [cf,idx] = sort(cf); cx=cx(idx,:);
            
            % End of Inner Loop for Competitive Evolution of Simplexes
        end;
        % Replace the complex back into the population;
        x(k2,:) = cx(k1,:);
        xf(k2) = cf(k1);
        
        % End of Loop on Complex Evolution;
    end;
    
    % Shuffled the complexes;
    [xf,idx] = sort(xf); x=x(idx,:);
    PX=x; PF=xf;
    
    % Record the best and worst points;
    bestx=x(1,:); bestf=xf(1);
    worstx=x(npt,:); worstf=xf(npt);
    BESTX=[BESTX;bestx]; BESTF=[BESTF;bestf];ICALL=[ICALL;icall];
    % Compute the standard deviation for each parameter
    xnstd=std(x);
    % Computes the normalized geometric range of the parameters
    gnrng=exp(mean(log((max(x)-min(x))./bound)));
    %     disp(['Evolution Loop: ' num2str(nloop) '  - Trial - ' num2str(icall)]);
    %     disp(['BESTF  : ' num2str(bestf)]);
    %     disp(['BESTX  : ' num2str(bestx)]);
    %     disp(['WORSTF : ' num2str(worstf)]);
    %     disp(['WORSTX : ' num2str(worstx)]);
    %     disp(' ');
    % Check for convergency;
    if debug
    if icall >= maxn;
        disp('*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE LIMIT');
        disp(['ON THE MAXIMUM NUMBER OF TRIALS ' num2str(maxn) ' HAS BEEN EXCEEDED!']);
    end;
    if gnrng < peps;
        disp('THE POPULATION HAS CONVERGED TO A PRESPECIFIED SMALL PARAMETER SPACE');
    end;
    end
    criter=[criter;bestf];
    if (nloop >= kstop);
        criter_change=abs(criter(nloop)-criter(nloop-kstop+1))*100;
        criter_change=criter_change/mean(abs(criter(nloop-kstop+1:nloop)));
        if criter_change < pcento;
            if debug
            disp(['THE BEST POINT HAS IMPROVED IN LAST ' num2str(kstop) ' LOOPS BY ', ...
                'LESS THAN THE THRESHOLD ' num2str(pcento) '%']);
            disp('CONVERGENCY HAS ACHIEVED BASED ON OBJECTIVE FUNCTION CRITERIA!!!')
            end
        end;
    end;
    
    % End of the Outer Loops
end;
if debug
disp(['SEARCH WAS STOPPED AT TRIAL NUMBER: ' num2str(icall)]);
disp(['NORMALIZED GEOMETRIC RANGE = ' num2str(gnrng)]);
disp(['THE BEST POINT HAS IMPROVED IN LAST ' num2str(kstop) ' LOOPS BY ', ...
    num2str(criter_change) '%']);
end

for p=1:par_size
    param_opt.(fn{select_par(p)}) = bestx(:,p);
end
% END of Subroutine sceua
return;
