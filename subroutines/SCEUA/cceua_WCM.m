%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [snew,fnew,icall]=cceua_WCM(s,sf,bl,bu,icall,maxn,fwd_in,obs,s_0,p_constraint,OF_param,param_0,select_par)

%  This is the subroutine for generating a new point in a simplex
%  Michel Bechtold 20 Jun 2021, adapted for use with WCM
%  Gabrielle De Lannoy 14 Mar 2023, updated calls to OF
%
%   s(.,.) = the sorted simplex in order of increasing function values
%   s(.) = function values in increasing order
%
% LIST OF LOCAL VARIABLES
%   sb(.) = the best point of the simplex
%   sw(.) = the worst point of the simplex
%   w2(.) = the second worst point of the simplex
%   fw = function value of the worst point
%   ce(.) = the centroid of the simplex excluding wo
%   snew(.) = new point generated from the simplex
%   iviol = flag indicating if constraints are violated
%         = 1 , yes
%         = 0 , no

fn          = fieldnames(param_0.m);
par_size    = length(select_par);       % number of parameters
param_0_tmp = param_0;                  % local control vector

[nps,nopt]=size(s);
n = nps;
m = nopt;
alpha = 1.0;
beta = 0.5;
% Assign the best and worst points:
sb=s(1,:); fb=sf(1);
sw=s(n,:); fw=sf(n);
% Compute the centroid of the simplex excluding the worst point:
ce=mean(s(1:n-1,:));
% Attempt a reflection point
snew = ce + alpha*(ce-sw);
% Check if is outside the bounds:
ibound=0;
s1=snew-bl; idx=find(s1<0); if ~isempty(idx); ibound=1; end;
s1=bu-snew; idx=find(s1<0); if ~isempty(idx); ibound=2; end;
if ibound >=1;
    snew = bl + rand(1,nopt).*(bu-bl);
end;
% Water Cloud Model
for p=1:par_size
    param_0_tmp.m.(fn{select_par(p)})  = snew(p);
end
if (strcmp(OF_param.OF_type,'SSE'))
    [OF_tmp, OF_sse_tmp, OF_p_tmp] = ...
        CalWCM_SSE(fwd_in, ...
        param_0_tmp.m, obs, p_constraint, ...
        OF_param,select_par);
else
    error(['not ready for OF-option ',OF_param.OF_type])
end    
fnew = OF_tmp;
% Reflection failed; now attempt a contraction point:
if fnew > fw;
    snew = sw + beta*(ce-sw);
    % Water Cloud Model
    for p=1:par_size
        param_0_tmp.m.(fn{select_par(p)})  = snew(p);
    end
    [OF_tmp, OF_sse_tmp, OF_p_tmp] = ...
        CalWCM_SSE(fwd_in, ...
        param_0_tmp.m, obs, p_constraint, ...
        OF_param,select_par);
    fnew = OF_tmp;
    % Both reflection and contraction have failed, attempt a random point;
    if fnew > fw;
        snew = bl + rand(1,nopt).*(bu-bl);
        % Water Cloud Model
        for p=1:par_size
            param_0_tmp.m.(fn{select_par(p)})  = snew(p);
        end
        [OF_tmp, OF_sse_tmp, OF_p_tmp] = ...
            CalWCM_SSE(fwd_in, ...
            param_0_tmp.m, obs, p_constraint, ...
            OF_param,select_par);
        fnew = OF_tmp;
    end;
end;
% END OF CCE
return;
