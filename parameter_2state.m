%% Evaluating residual for parameter estimation
% Calls fcn_2state.m

function y = parameter_2state(cdb1b2pb,lam1,lam2,N1,N2,b,K,b_eff,h,tmax,DT,...
                              data,pos_y_initial)

% Extract parameters from vector of estimated parameters
% This vector can be adjusted to reflect the parameters to be estimated
% Here an initial estimation that includes parameters p and bt is
% illustrated

c  = cdb1b2pb(1);    % speed
d  = cdb1b2pb(2);    % diffusion coefficient
b1 = cdb1b2pb(3);    % rate from movement to diffusion
b2 = cdb1b2pb(4);    % rate from diffusion to movement
p  = cdb1b2pb(5);    % fraction of particles starting in movement
bt = cdb1b2pb(6);    % parameter bt

%% Compute fluorescence recovery and residual with data

U    = fcn_2state(bt,c,d,b1,b2,p,lam1,lam2,N1,N2,b,K,b_eff,h,tmax,DT,...
                  pos_y_initial);
U    = U(1:end); 
y    = U-data;   % Need this setup for Multisearch with least squares norm

end
