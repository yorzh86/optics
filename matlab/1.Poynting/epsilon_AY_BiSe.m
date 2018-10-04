function [ep_cE] = epsilon_AY_BiSe(om)
    c0 = 299792458;
   
    wl = 2*pi*c0/om*1e6;
    w = 1.24/wl;
    
    ep_inf = 2.55;
    
    w0 = 1.0585527;
    wp = 2.63263661;
    gamma = 0.13258352;
    
    Re_ep_cE = ep_inf + wp^2*(w0^2-w^2)/((w0^2-w^2)^2+w^2*gamma^2);
    Im_ep_cE = wp^2*gamma*w/((w0^2-w^2)^2+w^2*gamma^2);
    
    ep_cE = Re_ep_cE + 1i*Im_ep_cE;

end
