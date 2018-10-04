function [sigma, epsilon] = BiSe_conductivity(om, T, mu)

h = 1.054571596e-34;           %Planck's constant over 2*pi [J-s]
hh = h*2*pi;
heV = 6.5821192815e-16;        %[eV-s]
kb = 1.380648813e-23;          %Boltzmann constant [J/K]
kbeV = 8.617332478e-5;         %Boltzmann constant [eV/K]
ee = 1.60217656535e-19;          %elementary electric charge [C] or [J/ev]
ep0 = 8.854187817e-12;          %electric permittivity
c0 = 2.99792458e+8;
om_Oph=3.0386e+14;      %Optical phonon frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATERIAL DEFINITIONS AND DOMAIN
wl = 2*pi*c0/om*1e6;

thickness = 0.92e-9;
tau = 100e-15;  %Qu et al 2010 Science

sigma_Drude = 1i/(om+1i/tau)*ee^2/(pi*h^2)*abs(ee*mu);

xi = (0:0.001:1e2)*mu;
G_hom2 = sinh(heV*om/2/(kbeV*T))/(cosh(mu/(kbeV*T))+cosh(heV*om/2/(kbeV*T)));

for ii=1:length(xi)
    G(ii) = sinh(xi(ii)/(kbeV*T))/(cosh(mu/(kbeV*T))+cosh(xi(ii)/(kbeV*T)));
        if isnan(G(ii))
            G(ii) = 1;
        else
            G(ii) = G(ii);
        end
    Int_term(ii) = (G(ii)-G_hom2)/((heV*om)^2-4*xi(ii)^2);
end

% Integration = (max(xi)-min(xi))/(3*(length(xi)-1))*(Int_term(1)+Int_term(end)+4*sum(Int_term(2:2:(end-1)))+2*sum(Int_term(3:2:(end-2))));
Integration = trapz(xi,Int_term);
Second_Inter = 4*heV*om/pi*Integration;
sigma_Inter = ee^2/(4*h)*(G_hom2+1i*Second_Inter);
    
sigma = sigma_Drude + sigma_Inter;

epsilon=1+1i*sigma/(om*ep0*thickness); %important 1+1i sigma .... 
end
