function epsilon_Ag=epsilon_Ag(om)
c0=2.99792458e+8;       %m/s
wl=2*pi*c0/om*1e6;     %um
me = 9.10938291e-31;
ee = 1.60217657e-19;
ep0 = 8.85418782e-12;

% e_inf=4.08598;
% w_p=1.3316e16   %rad/s
% gama=1.1308e14   %rad/s
% eps_Ag = e_inf-w_p^2./(w.*(w+i*gama)); 

m0 = 0.96;
N = 5.85e28;
om_p = N*ee^2/(me*m0*ep0);
gamma = 1/31e-15;

gamma = 2.73e13;
om_p = 1.39e16;
epsilon_Ag = 1-om_p^2/(om^2+1i*gamma*om);
return

    
     
        
        
        
        
        
        
        
        
        
        
        
        
    