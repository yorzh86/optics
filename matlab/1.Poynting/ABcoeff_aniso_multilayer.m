% Richard Z. Zhang 05/08/2015
% This is for any # layers of anisotropic permittivity slabs.
% NOT FOR TILTED ANISO MEDIA
% ONLY TM WAVE
% Returns A and B coefficients for Poynting vector calculations
% using another function. Based on derivations from BJ Lee (2006)

function [A,B]=ABcoeff_aniso_multilayer(diff, layers, wl, ep_O, ep_E, kx)
c0 = 2.99792458e8;
ep0 = 8.854187817e-12;
om = 2*pi*c0/wl*1.0e6;
k0 = om/c0;

for ii=1:layers
    kz_p(ii) = sqrt(k0^2*ep_O(ii)-kx^2*ep_O(ii)/ep_E(ii));
    kz_n(ii) = -kz_p(ii);
    if (imag(kz_p(ii))<0)
        temp=kz_p(ii);
        kz_p(ii)=kz_n(ii);
        kz_n(ii)=temp;
    end
    z(ii) = 1/(om*ep0)*kz_p(ii)/ep_O(ii);
end

M = cell(1,layers-1);
M_prod = cell(1,layers);
M_prod{layers} = eye(2);
for jj=(layers-1):-1:1
    r_ij = (z(jj)-z(jj+1))/(z(jj)+z(jj+1));
    r_ji = (z(jj+1)-z(jj))/(z(jj)+z(jj+1));
    t_ij = (2*z(jj))/(z(jj)+z(jj+1));
    t_ji = (2*z(jj+1))/(z(jj)+z(jj+1));
    DD = 1/t_ij*[1, -r_ji; r_ij, (t_ij*t_ji-r_ij*r_ji)];
    P = [exp(-1i*kz_p(jj)*diff(jj)), 0.0; 0.0, exp(-1i*kz_n(jj)*diff(jj))];
    M{jj} = P*DD;
    M_prod{jj} = M{jj}*M_prod{jj+1};
end

%The transfer matrix for p-wave (for the Magnetic field)
A = zeros(1,layers);
B = zeros(1,layers);
B(layers) = 0;
A(layers) = 1/M_prod{1}(1,1);
B(1) = M_prod{1}(2,1)*A(layers);
A(1) = 1;
for ab=(layers-1):-1:1
    A(ab) = M_prod{ab}(1,1)*A(layers);
    B(ab) = M_prod{ab}(2,1)*A(layers);
end
return

