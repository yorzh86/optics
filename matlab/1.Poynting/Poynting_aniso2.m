% Richard Z. Zhang 4/22/2015
% To obtain the position-dependent Poynting vector components:
% S_x and S_z for both polarizations
% For any number of layers.

function [S_xTE, S_zTE, S_xTM, S_zTM] = Poynting_aniso2(A_TE, B_TE, A_TM, B_TM, z_norm, diff, diff_norm, k_x, ep_o, ep_e, wl)
mu0 = 4*pi/(10^7);
epsi0 = 8.854187817e-12;
c0 = 2.99792458e8;

om = 2*pi*c0/(wl/1.0e6);
k0 = om/c0;

z = z_norm*sum(diff);
z_0 = [0, sum(diff)*diff_norm];

for ii = 1:length(diff_norm)-1
    if (z_norm < diff_norm(1))
        ind = 1;
    elseif (z_norm < diff_norm(ii+1))&&(z_norm >= diff_norm(ii))
        ind = ii+1;
    elseif (z_norm >= diff_norm(end))
        ind = length(diff_norm)+1;
    end
end

kzTM(ind)=sqrt(k0^2*ep_o(ind)-ep_o(ind)/ep_e(ind)*k_x^2);
if (imag(kzTM(ind))<0)
    kzTM(ind)=-sqrt(k0^2*ep_o(ind)-ep_o(ind)/ep_e(ind)*k_x^2);
end

kzTE(ind) = sqrt(ep_o(ind)*k0^2-k_x^2);
if (imag(kzTE(ind))<0)
    kzTE(ind) = -sqrt(ep_o(ind)*k0^2-k_x^2);
end

% For TM waves only:
Hy = (A_TM(ind)*exp(1i*kzTM(ind)*(z-z_0(ind)))+B_TM(ind)*exp(-1i*kzTM(ind)*(z-z_0(ind))))*exp(1i*k_x);
Hy_conj = conj(Hy);
dHy_dz = (kzTM(ind)*A_TM(ind)*exp(1i*kzTM(ind)*(z-z_0(ind)))-kzTM(ind)*B_TM(ind)*exp(-1i*kzTM(ind)*(z-z_0(ind))))*exp(1i*k_x);
dHy_dx = k_x*Hy;
Ex = 1/(om*epsi0*ep_o(ind)*ep_e(ind))*(ep_e(ind)*dHy_dz);
Ez = -1/(om*epsi0*ep_o(ind)*ep_e(ind))*(ep_o(ind)*dHy_dx);
S_xTM = -0.5*real(Ez*Hy_conj);
S_zTM = 0.5*real(Ex*Hy_conj);

% For TE waves only:
Ey = (A_TE(ind)*exp(1i*kzTE(ind)*(z-z_0(ind)))+B_TE(ind)*exp(-1i*kzTE(ind)*(z-z_0(ind))))*exp(1i*k_x);
dEy_dz = (1i*kzTE(ind)*A_TE(ind)*exp(1i*kzTE(ind)*(z-z_0(ind)))-1i*kzTE(ind)*B_TE(ind)*exp(-1i*kzTE(ind)*(z-z_0(ind))))*exp(1i*k_x);
dEy_dx = 1i*k_x*(A_TE(ind)*exp(1i*kzTE(ind)*(z-z_0(ind)))+B_TE(ind)*exp(-1i*kzTE(ind)*(z-z_0(ind))))*exp(1i*k_x);
Hx = -1/(1i*om*mu0)*dEy_dz;
Hz = 1/(1i*om*mu0)*dEy_dx;
Hx_conj = conj(Hx);
Hz_conj = conj(Hz);
S_xTE = 0.5*real(Ey*Hz_conj);
S_zTE = -0.5*real(Ey*Hx_conj);
return
