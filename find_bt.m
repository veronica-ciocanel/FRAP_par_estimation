%% Function that finds parameter bt in the Exponential of Gaussian 
% initial condition to ensure that FRAP data starts with the initial 
% fluorescence value a

function bt  = find_bt(K,b_eff,b,lam1,lam2,N1,N2,pos_y_initial,a)

% Spatial grid
x = (lam1/N1)*(1:N1)'; y = (lam2/N2)*(1:N2)'; [X,Y] = ndgrid(x,y);           
radius = sqrt((X-lam1/2).^2+(Y-pos_y_initial).^2);

% Initial condition given by Exp(Gaussian), approximates FRAP bleach spot
uv = exp(-K.*exp(-2.*radius.^2/b_eff^2)); 
uv_bleach1 = (radius<=b).*real(uv);               % u+v over bleach spot, 0 elsewhere
bt = a/(trapz(x,trapz(y,uv_bleach1,2),1)/(pi*b^2));

end


