% Epsilon of Bi2Se3 from
% M.S. Wolf, Infrared and Optical Studies of Topological Insulators Bi2Se3, Bi2Te3 and Sb2Te3
% R Zhang 110217

function [Re_ep, Im_ep] = epsilon_BiSe(om)
    c0 = 299792458;
   
    wl = 2*pi*c0/om*1e6;
    w = 1E4/wl;
    
    ep_inf = 1;
    
    % Drude (cm^-1)
    w_D = 908.66;
    gamma_D = 7.43;
    
    ep_Drude = -w_D^2/(w^2+1i*w*gamma_D);
    
    %Alpha Lorentz (cm^-1)
    w0_A = 63.03;
    wp_A = 675.9;
    gamma_A = 17.5;
    
    ep_Alpha = wp_A^2/(w0_A^2-w^2-1i*w*gamma_A);
    
    %Beta Lorentz (cm^-1)
    w0_B = 126.94;
    wp_B = 100;
    gamma_B = 10;
    
    ep_Beta = wp_B^2/(w0_B^2-w^2-1i*w*gamma_B);
    
    %Omega Lorentz (cm^-1)
    w0_g = 2029.5;
    wp_g = 11249;
    gamma_g = 3920.5;
    
    ep_gap = wp_g^2/(w0_g^2-w^2-1i*w*gamma_g);
    
    ep = ep_inf + ep_Drude + ep_Alpha + ep_Beta + ep_gap;
    Re_ep = real(ep);
    Im_ep = imag(ep);
return;