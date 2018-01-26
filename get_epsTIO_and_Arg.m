wl = [0.4, 0.5, 0.6, 0.7, 0.8];


c0 = 2.99792458e+8;            %speed of light in vacuum
ep0 = 8.854187817e-12;

for ii = 1:length(wl)
    om = 2*pi*c0./wl(ii)*1.0e6;
    [ep_TiO, ep_TiE] = epsilon_TiO2(om);
    X = ["wavelength:", wl(ii), 'eps_TiO2_E:', ep_TiE];
    Y = ["wavelength:", wl(ii), 'eps_TiO2_O:', ep_TiO];
    disp(X);
    disp(Y);
end

c0=2.99792458e+8;       %m/s
%wl=2*pi*c0/om*1e6;     %um
me = 9.10938291e-31;
ee = 1.60217657e-19;
ep0 = 8.85418782e-12;

m0 = 0.96;
N = 5.85e28;
om_p = N*ee^2/(me*m0*ep0);
gamma = 1/31e-15;

gamma = 2.73e13;
om_p = 1.39e16;
for ii = 1:length(wl)
    om = 2*pi*c0./wl(ii)*1.0e6;
    epsilon_Ag = 1-om_p^2/(om^2+1i*gamma*om);
    X = ["wavelength:", wl(ii), 'eps_Ag:', epsilon_Ag];
    disp(X)
end
