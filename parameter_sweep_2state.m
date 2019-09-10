%% Parameter sweep for the 2-state model 
% Helps choose an initial condition for the parameter estimation 
% simulations: here, we use the movement-diffusion 2-state system for 
% WT data for an average FRAP dataset
% Uses 'average_dataset.mat'
% Calls fcn_2state.m and find_bt.m
% Main code for parameter sweeps

% c: speed of particles that are actively transported
% D: diffusion coefficient of particles that are free to diffuse
% b1: reaction rate from movement state to diffusion state
% b2: reaction rate from diffusion state to movement state

close all;
clear all; 

type    = 'WT'; % type of oocyte considered
oocyte_names = [222,223,234,240,241];
region  = 1;  % region
use_par = 0;  % index is 1 if using parallel computing

% Initial condition parameters - from fitting postbleach intensity profile
K      = 0.88; % average
b_eff  = 10.9; % average

% Spatial grid and other parameters
lam1 = 40; lam2 = 60; N1 = 64; N2 = 64; % domain/grid size
b   = 2.5;           % radius of FRAP bleach spot
h   = 0.1;           % timestep 
pos_y_initial = lam2/2; % position on the y-axis where the bleach is 
                        % initialized

% Time parameters
tmax     = 200;     % target time  
DT       = 5;       % time increment
time_vect = [0:DT:tmax]; % time vector

% Load average FRAP dataset
a_temp = load('average_dataset.mat'); 
data = a_temp.data;    
a    = data(1);

% Find parameter bt which ensures the intensity in the initial bleach spot 
% matches the FRAP data
bt  = find_bt(K,b_eff,b,lam1,lam2,N1,N2,pos_y_initial,a);

dir_name = ['FRAP_2state'];
test     = exist(dir_name,'dir');
if test(1,1) == 0
	mkdir(dir_name);
end
 
% Adapt this parameter sweep for c as necessary
c_col   = [0.0001,0.001,0.01:0.01:0.09,0.1:0.1:0.5];  % (16 values) 

NT = size(c_col,2);  % number of trials/values of speed c  

% If using parallel computing, start parallel pool, with one cluster
% for each value of speed c tested
if use_par == 1
    myCluster = parcluster('local');
    delete(myCluster.Jobs);
    poolobj = parpool(NT);  
end

for tr = 1:NT   % Use for if parallel computing is not available
%parfor tr = 1:NT  % Use parfor when using parallel computing

c = c_col(tr);

tic;    

% Adapt this parameter sweep for D, b1 and b2 as necessary
% parameters tested   
D_log1  = [0.005,0.01,0.05,0.1,0.5,1];  
b1_log  = -6:1:2;  
b1_log1 = 10.^b1_log;
b2_log  = -6:1:2;  
b2_log1 = 10.^b2_log;
p_123   = [0,0.25,0.5,0.75,1];

i = 1; % index
nr_rows   = length(D_log1)*length(b1_log1)*length(b2_log1)*length(p_123);
table_val = zeros(nr_rows,7); % matrix that stores the results

% Evaluate the main function at all the parameter combinations and store
% the results in table
    for D = D_log1
        for b1 = b1_log1
            for b2 = b2_log1
                for p = p_123

                    [flor] = fcn_2state(bt,c,D,b1,b2,p,lam1,lam2,N1,N2,b,...
                                        K,b_eff,h,tmax,DT,pos_y_initial);
                    % residual between data generated with given set of
                    % parameters and actual FRAP data
                    yval   = sqrt(sum((flor - data).^2));   
                    
                    % store parameters and residual in table
                    table_val(i,:) = [c D b1 b2 p bt yval];
                    i = i+1;
                    
                end
            end
        end
    end

elapsed_time  = toc; 
    
% Save results in as many files as parameters of speed c tested
file_name_new = [dir_name '/flor_2state_' type '_ROI' num2str(region)...
                '_c' num2str(c) '.mat'];
parsave1(file_name_new,table_val,flor,type,region,D_log1,b1_log1,b2_log1,...
         p_123,elapsed_time);
   
end

if use_par == 1
    delete(poolobj);
end
who
if use_par == 1
    clear myCluster;
end
