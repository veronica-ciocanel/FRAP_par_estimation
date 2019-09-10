%% Function that generated the FRAP curve for given parameters
% Example run: 
% flor =
% fcn_2state(30,0.12,1.63,0.0005,0.02,0.07,40,60,64,64,2.5,0.88,10.9,...
%                  0.1,200,5,50);
% figure(1)
% plot(0:5:200,flor,'*','Linewidth',2)

function flor = fcn_2state(bt,c,D,b1,b2,p,lam1,lam2,N1,N2,b,K,b_eff,h,tmax,...
                           DT,pos_y_initial)

% bt: this parameter controls the magnitude of the initial conditions so
% that the intensity in the bleach spot matches the FRAP data
% c: speed of particles in active transport/movement state
% D: diffusion coefficient of particles in diffusion state
% b1: reaction rate from movement state to diffusion state
% b2: reaction rate from diffusion state to movement state
% p: proportion of particles undergoing active transport initially
% lam1: length of horizontal (x) dimension of the domain considered
% lam2: length of vertical (y) dimension of the domain considered
% N1: number of grid points in x direction
% N2: number of grid points in y direction
% b: radius of FRAP bleach spot
% K: bleaching depth of bleach spot from fit of postbleach profile
% b_eff: effective radius of bleach spot from fit of postbleach profile
% h: time step for time integration
% tmax: time length of the FRAP experiment
% DT: interval of time between FRAP readings
% pos_y_initial: position on the y-axis where the bleach is initialized

% The code below is adapted from "Solving reaction-diffusion 
% equations 10 times faster" by Kassam and "Fourth-order time-stepping
% for stiff PDE" by Kassam and Trefethen

%========================= GENERIC SET UP ========================= 
% Spatial grid
x = (lam1/N1)*(1:N1)'; y = (lam2/N2)*(1:N2)'; [X,Y] = ndgrid(x,y);           
radius = sqrt((X-lam1/2).^2+(Y-pos_y_initial).^2);

% Initial conditions given by an exponential of a Gaussian
% Approximates FRAP bleach spot
u = p*(bt.*exp(-K.*exp(-2.*radius.^2/b_eff^2))); % concentration of...
% particles in active transport/movement state    
u_ft = fftn(u); 
v = (1-p)*(bt.*exp(-K.*exp(-2.*radius.^2/b_eff^2))); % concentration of...
% particles in diffusion state
v_ft = fftn(v); 

ntimes    = round(tmax/DT);
flor      = zeros(1,ntimes);            % fluorescence recovery vector
uv_bleach1= (radius<=b).*real(u+v);     % u+v over bleach spot, 0 elsewhere
% intensity over bleach spot:
flor(1)   = trapz(x,trapz(y,uv_bleach1,2),1)/(pi*b^2);

% Precompute various ETDRK4 matrix quantities
k1        = [0:N1/2-1 0 -N1/2+1:-1]'/(lam1/(2*pi)); % wave numbers
k2        = [0:N2/2-1 0 -N2/2+1:-1]'/(lam2/(2*pi)); % wave numbers
[xi,eta]  = ndgrid(k1,k2);             % 2D wave numbers. 
L_u       = c*1i*eta - b1;             % movement/convection down (y-axis)  
L_v       = -D*(eta.^2+xi.^2) - b2 ;   % 2D Laplacian + reaction term.

Fr1       = false(N1,1);           %High frequencies for de-aliasing
Fr1(N1/2+1-round(N1/6) : N1/2+round(N1/6))=1;   
Fr2       = false(N2,1);           %High frequencies for de-aliasing
Fr2(N2/2+1-round(N2/6) : N2/2+round(N2/6))=1; 
[alxi,aleta] = ndgrid(Fr1,Fr2);  
ind          = alxi | aleta; 

%=============== PRECOMPUTING ETDRK4 COEFFS ===================== 
E_u = exp(h*L_u); E2_u = exp(h*L_u/2);
E_v = exp(h*L_v); E2_v = exp(h*L_v/2); 
M   = 16; % number of points for complex mean 
r   = exp(1i*pi*((1:M)-0.5)/M); % roots of unity 
L_u = L_u(:); LR_u = h*L_u(:,ones(M,1))+r(ones(N1*N2,1),:);  
L_v = L_v(:); LR_v = h*L_v(:,ones(M,1))+r(ones(N1*N2,1),:);

Q_u = h*real(mean( (exp(LR_u/2)-1)./LR_u ,2));
Q_v = h*real(mean( (exp(LR_v/2)-1)./LR_v ,2));

f1_u = h*real(mean( (-4-LR_u+exp(LR_u).*(4-3*LR_u+LR_u.^2))./LR_u.^3 ,2));
f1_v = h*real(mean( (-4-LR_v+exp(LR_v).*(4-3*LR_v+LR_v.^2))./LR_v.^3 ,2));

f2_u = h*real(mean( (4+2*LR_u+exp(LR_u).*(-4+2*LR_u))./LR_u.^3 ,2)); 
f2_v = h*real(mean( (4+2*LR_v+exp(LR_v).*(-4+2*LR_v))./LR_v.^3 ,2)); 

f3_u = h*real(mean( (-4-3*LR_u-LR_u.^2+exp(LR_u).*(4-LR_u))./LR_u.^3 ,2)); 
f3_v = h*real(mean( (-4-3*LR_v-LR_v.^2+exp(LR_v).*(4-LR_v))./LR_v.^3 ,2)); 

f1_u = reshape(f1_u,N1,N2); f2_u=reshape(f2_u,N1,N2); f3_u=reshape(f3_u,N1,N2);
f1_v = reshape(f1_v,N1,N2); f2_v=reshape(f2_v,N1,N2); f3_v=reshape(f3_v,N1,N2);

Q_u = reshape(Q_u,N1,N2); clear LR_u 
Q_v = reshape(Q_v,N1,N2); clear LR_v


%==================== TIME STEPPING LOOP ======================= 
nmax   = round(tmax/h); 
i      = 2;

for n = 1:nmax 
    t  = n*h; %**Nonlinear terms are evaluated in physical space** 
    Nu = fftn( g_u(ifftn(v_ft),b2) ); %Nonlinear evaluation. g(u,*)
    Nv = fftn( g_v(ifftn(u_ft),b1) ); 
    
    a_u = E2_u.*u_ft + Q_u.*Nu; %Coefficient a in ETDRK formula 
    a_v = E2_v.*v_ft + Q_v.*Nv; 
    
    Na_u = fftn( g_u(ifftn(a_v),b2) ); %Nonlinear evaluation. g(a,*)
    Na_v = fftn( g_v(ifftn(a_u),b1) );
    
    b_u = E2_u.*u_ft + Q_u.*Na_u; %Coefficient b in ETDRK formula 
    b_v = E2_v.*v_ft + Q_v.*Na_v;
    
    Nb_u = fftn( g_u(ifftn(b_v),b2) ); %Nonlinear evaluation. g(b,*) 
    Nb_v = fftn( g_v(ifftn(b_u),b1) );
    
    c_u = E2_u.*a_u + Q_u.*(2*Nb_u-Nu); %Coefficient c in ETDRK formula 
    c_v = E2_v.*a_v + Q_v.*(2*Nb_v-Nv);
    
    Nc_u = fftn( g_u(ifftn(c_v),b2) ); %Nonlinear evaluation. g(c,*)
    Nc_v = fftn( g_v(ifftn(c_u),b1) );
    
    u_ft = E_u.*u_ft + Nu.*f1_u + (Na_u+Nb_u).*f2_u + Nc_u.*f3_u; %update
    v_ft = E_v.*v_ft + Nv.*f1_v + (Na_v+Nb_v).*f2_v + Nc_v.*f3_v;
    
    u_ft(ind) = 0; % High frequency removal --- de-aliasing
    v_ft(ind) = 0;
      
    if mod(t,DT)==0 
        u = ifftn(u_ft); 
        v = ifftn(v_ft);
        
        %=============== FRAP DATA ===============
        uv_bleach = (radius<=b).*real(u+v);% u+v over bleach spot, 0 elsewhere
        % intensity over bleach spot
        flor(i)   = trapz(x,trapz(y,uv_bleach,2),1)/(pi*b^2); 
        i         = i+1;                   % update index
       
        % It may also be useful to add a plot of the solutions u, v, and u+v 
        % as a function of time here to ensure that the domain chosen is
        % large enough and that discretization parameters are fine enough
        % to avoid numerical instabilities. This is also a useful way to
        % run the PDEs with estimated parameters and check results.
    end 
    
end
        
end


function y = g_u(v,b2)
y = b2*v;  % function for u
end

function y = g_v(u,b1)
y = b1*u;  % function for v
end

