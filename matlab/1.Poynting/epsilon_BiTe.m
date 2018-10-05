% Epsilon of Bi2Se3 from
% M.S. Wolf, Infrared and Optical Studies of Topological Insulators Bi2Se3, Bi2Te3 and Sb2Te3
% A. Yorzh 1052018

function [Re_ep, Im_ep] = epsilon_BiTe(om)
    c0 = 299792458;
   
    wl = 2*pi*c0/om*1e6;
    w = 1E4/wl;
    
    ep_inf = 1;
    
    % Drude (cm^-1)
    w_D = 5651.5;
    gamma_D = 111.86;
    
    ep_Drude = -w_D^2/(w^2+1i*w*gamma_D);
    
    %Alpha Lorentz (cm^-1)
    w0_A = 8386.6;
    wp_A = 66024;
    gamma_A = 10260;
    
    ep_Alpha = wp_A^2/(w0_A^2-w^2-1i*w*gamma_A);
    
    ep = ep_inf + ep_Drude + ep_Alpha;
    Re_ep = real(ep);
    Im_ep = imag(ep);
return;