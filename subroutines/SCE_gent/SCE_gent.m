function [param_opt] = ... 
       SCE_gent(param_0, select_par, fwd_in, obs)

% Gabrielle De Lannoy - 13 March 2023: function extracted
%      Code given to Shannon, used by Lievens et al. (WCM calib paper)
%==========================================================================
   
   fn          = fieldnames(param_0.m);
   par_size    = length(select_par);       

   X0 = []; LB = []; UB = [];
   for p=1:par_size
        X0(p) = param_0.m.(fn{select_par(p)});
        LB(p) = param_0.r.(fn{select_par(p)})(1);
        UB(p) = param_0.r.(fn{select_par(p)})(2);
   end 
   mat=[fwd_in.sm fwd_in.veg obs]; 
   mat_cell = mat2cell(mat',[3]);
   
   %CALWCM.m is hardwired to optimize all 4 parameters, 
   %assuming theta = 0 deg, and
   %using an ad hoc combination of OF = RMSE + 0.1*sum(prior dev)
   its = 2; %questionable, b/c not actually used to optimize
   ns0s = 1; 
   for it=1:its
    for s=1:ns0s
             [X,~,~,~] = SCE('CALWCM',X0',LB',UB',[ ],mat_cell(s));
             display(['SCE: ',num2str(it),' ',num2str(s)])
             for p=1:par_size
                param_opt_tmp.(fn{select_par(p)})(s,it) = X(p);
                display(['SCE par: ',num2str(X(p))])
             end 
    end
   end
   for p=1:par_size
       param_opt.(fn{select_par(p)}) = param_opt_tmp.(fn{select_par(p)})(1,1);
   end 

end
%==========================================EOF=============================
