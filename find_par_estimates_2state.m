%% Find parameters in the advection-diffusion 2-state model PDE using Multisearch 
% Uses 'average_dataset.mat'
% Calls parameter_2state.m, fcn_2state.m, and find_bt.m
% Main code to run for parameter estimation

close all;
clear all; 

type    = 'WT'; % type of oocyte considered
region  = 1;    % region
use_par = 0;    % 0 if no parallel processing, 1 for parallel processing
rand_st = 0;    % 1 if using random start points, 0 if custom start points
NT      = 5;    % number of cores used if run in parallel 
k_pcts  = 5;    % number of start points (only relevant for random start points)

dir_name = ['Multisearch_2state'];
test     = exist(dir_name,'dir');
if test(1,1) == 0
	mkdir(dir_name);
end

% Initial condition parameters - from fitting postbleach intensity profile
K      = 0.88; 
b_eff  = 10.9; 

% Initial guess: needed for MultiSearch
c_star  = 0.03;
D_star  = 0.5;
b1_star = 0.005;
b2_star = 0.005;
p_star  = 0.25;
h       = 0.1;  % timestep

tic;    

%% Time parameters
tmax     = 200;     % target time 
DT       = 5;       % time increment
time_vect = [0:DT:tmax]; % time vector

%% TolX, TolFun (tolerances for lsqnonlin optimization)   
tol_nr1 = 1e-6;  
tol_nr2 = 1e-6;

%% Maximum evals/iters
max_fun_evals = 1e+10;
max_iters     = 1e+10;
 
%% Parallel Processing
if use_par == 1          
    myCluster = parcluster('local');
    delete(myCluster.Jobs);  
    parpool(NT);
end

%% Grid and other constants
lam1 = 40; lam2 = 60; N1 = 64; N2 = 64; % domain/grid size
b    = 2.5;          % original radius of FRAP bleach spot 
pos_y_initial = lam2/2; % position on the y-axis where bleach is initialized
lim_value = 200;     % value to penalize for unrealistic parameters

%% FRAP Data   
% Load average FRAP dataset
a_temp = load('average_dataset.mat'); 
data = a_temp.data;    
a    = data(1);

% Find parameter bt which ensures the intensity in the bleach spot 
% matches the initial postbleach intensity
bt      = find_bt(K,b_eff,b,lam1,lam2,N1,N2,pos_y_initial,a);
bt_star = bt;

%% Initial guess
initial_guess = [c_star D_star b1_star b2_star p_star bt_star];

%% MultiStart set-up
options = optimoptions(@lsqnonlin,'FinDiffType','central',...
          'MaxIter',max_iters,'MaxFunEvals',max_fun_evals);
% The bounds below may be made tighter to ensure convergence or expanded if
% estimates run out of bounds
LB = [0 0 0   0   0 0     ];  % lower bounds for parameters
UB = [1 5 10  10  1 bt+50 ];  % upper bounds for parameters

problem = createOptimProblem('lsqnonlin','objective',...
 @(cdb1b2pb)parameter_2state(cdb1b2pb,lam1,lam2,N1,N2,b,K,b_eff,h,tmax,DT,...
 data,pos_y_initial),'x0',initial_guess,'lb',LB,'ub',UB,'options',options);

if rand_st == 1
    % Random start points
    startpts = RandomStartPointSet('NumStartPoints',k_pcts);
else
    % Custom start points from ranking parameter sweep results in terms of
    % the smallest residual (start points come from parameter sweeps for 
    % different oocytes/trials as well as sweeps for average data)
    ptmatrix = [
                0.2      1       0.01     0.01     0       bt_star; 
                0.03     1       0.00001  0.01     0.5     bt_star; 
                0.01     0.5     0.001    0.01     0.25    bt_star;
                
                % average
                0.02     1       0.000001 0.001    0.25    bt_star;
                0.5      0.05    0.1      0.001    0.75    bt_star;
              
              ];  % each line of this matrix is an initial guess (8)
    startpts = CustomStartPointSet(ptmatrix);
end

if use_par==1
    ms = MultiStart('UseParallel',true,'StartPointsToRun','bounds',...
    'TolX',tol_nr1,'TolFun',tol_nr2);
else
    ms = MultiStart('UseParallel',false,'StartPointsToRun','bounds',...
    'TolX',tol_nr1,'TolFun',tol_nr2);
end

%% Run parameter estimation
[params_new,fval,exitflag,output,solutions] = run(ms,problem,startpts);
message = ['The optimized value of c, d, b1, b2 is ', num2str(params_new)];

elapsed_time = toc;

%% Save data
file_name_new = [dir_name '/data_ms_2state_' type '_ROI' num2str(region) '.mat'];
parsave1(file_name_new,params_new,fval,exitflag,output,solutions,elapsed_time,...
    message,a,K,b_eff,bt_star,p_star,c_star,D_star,b1_star,b2_star,...
    initial_guess,lim_value,options,tol_nr1,tol_nr2,region,type,h,...
    max_fun_evals,max_iters,time_vect,data,LB,UB,ms,problem,startpts,bt_star);

if use_par==1
    clear myCluster;
end

who

