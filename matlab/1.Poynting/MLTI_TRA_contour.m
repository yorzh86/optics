% Calculation of radiative coefficients of multilayer, and Poynting vector
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
kaishi=clock;
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0 = 2.99792458e+8;            %speed of light in vacuum
ep0 = 8.854187817e-12;
mu0 = 4*pi*1e-7;

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
wl = [0.5]; %1.5, 2.5
%wl = logspace(log10(0.4),log10(20),50); %need
%wl = linspace(0.4,20,50);
nu = 1e4./wl;
om = 2*pi*c0./wl.*1e6;

k0 = om./c0;
d_theta = 0.5;
%theta_i = 0:d_theta:90; %degrees %need
theta_i = [10, 30, 50, 70];

d_c = 0.92e-9;          %BiSe charge thickness
% d_c = 1.90e-9;        %BiTe charge thickness
d_d = 12e-9;
d_b = 10e-9-2*d_c;
per = 91;               %period add manually!!!

diff = zeros(1,4*per+2);
diff(2:4:end) = d_d;
diff(3:2:end-1) = d_c;
diff(4:4:end) = d_b;
diff(1) = 0;

step = 100;

z_min = -0.5;
z_max = 2;
slab = sum(diff);
z_range{1} = linspace(z_min,-1/step,step);
for zz = 2:length(diff)
    z_range{zz} = linspace(sum(diff(1:(zz-1)))/slab,sum(diff(1:zz))/slab-(diff(zz)/slab)/step,step);
end
z_range{length(diff)+1} = linspace(1,z_max,step);
z_norm = [z_range{1:end}];

origin = find(z_norm==0,1,'first');
backside = find(z_norm==1,1,'first');

stepEMT = 10;
diff_norm = cumsum(diff)./slab;
diff_EMT = [0, slab];
z_normEMT = [linspace(z_min, -1/stepEMT, stepEMT), linspace(0, 1-1/(stepEMT), stepEMT), linspace(1, z_max, stepEMT)];
dif_n = [0, 1]

slab = sum(diff);
diff_norm = cumsum(diff)./slab;
diff_EMT = [0, slab];

origEMT = find(z_normEMT==0,1,'first');

ff = d_b/(d_b+d_d);
f_d = d_d/(d_b+d_d+2*d_c);
f_b = d_b/(d_b+d_d+2*d_c);

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
T = 300;
mu = 0.189;
% mu = 0.194;

    %------------------------------------------
% poolobj = gcp('nocreate');  
% if isempty(poolobj)
%     parpool('local');
%     disp('parallel environment open successfully');
% else
%     disp('parallel pool already open. computing'); 
% end

for ii = 1:length(wl)
            %-DIELECTRIC EPSILON-----------------------------------
%         ep_realZnSe(ii) = real(4+1.9*wl(ii)^2/(wl(ii)^2-0.113));
        ep_ZnSe(ii) = 2.5^2;
        ep_ZnTe(ii) = 2.8^2;
        
        ep_dO(ii) = ep_ZnSe(ii);
        ep_dE(ii) = ep_ZnSe(ii);

        [Re_ep_TI, Im_ep_TI] = epsilon_BiSe(om(ii));
%         [Re_ep_TI, Im_ep_TI] = epsilon_BiTe(om(ii));
        ep_bO(ii) = Re_ep_TI+1i*Im_ep_TI;
        ep_bE(ii) = Re_ep_TI+1i*Im_ep_TI;

        %-METAL EPSILON-----------------------------------
        [sigma(ii),~]=BiSe_conductivity(om(ii),T,mu);
%         [sigma(ii),~]=BiTe_conductivity(om(ii),T,mu);
        %[sigma(ii),~]=Graphene_conductivity(om(ii),T,mu);

        [ep_Ec(ii)] = epsilon_AY_BiSe(om(ii));  
%         [ep_Ec(ii)] = epsilon_AY_BiTe(om(ii));  

        ep_c(ii)=1+1i*sigma(ii)/(om(ii)*ep0*d_c); %changed to 1+ 1i...

        ep_cO(ii) = ep_c(ii)+ep_bO(ii);
        ep_cE(ii) = ep_Ec(ii);

        ep_pO = [ep_dO(ii), ep_cO(ii), ep_bO(ii), ep_cO(ii)];    %period dielectric function
        ep_pE = [ep_dE(ii), ep_cE(ii), ep_bE(ii), ep_cE(ii)];    %these are assigned to diff

        ep1 = 1.0;      %vacuum in medium 1 and 3
        ep3 = 1.0;

        %------------------------------------------

        ep_MLO = repmat(ep_pO,1,per);
        ep_MLE = repmat(ep_pE,1,per);

        ep_allO = [ep1,ep_MLO,ep_dO(ii),ep3];
        ep_allE = [ep1,ep_MLE,ep_dE(ii),ep3];

        layers = length(ep_allO);

        ep_EMTO(ii) = ff*ep_bO(ii)+(1-ff)*ep_dO(ii);
        ep_EMTE(ii) = ep_bE(ii)*ep_dE(ii)/(ff*ep_dE(ii)+(1-ff)*ep_bE(ii));
        ep_iEMTO(ii) = f_b*ep_bO(ii)+f_d*ep_dO(ii)+(1-f_d-f_b)*ep_cO(ii);
        ep_iEMTE(ii) = 1/(f_b/ep_bE(ii)+f_d/ep_dE(ii)+(1-f_d-f_b)/ep_cE(ii));

        ep_3layerO = [ep1, ep_iEMTO(ii), ep3];
        ep_3layerE = [ep1, ep_iEMTE(ii), ep3];
        
    for ind = 1:length(theta_i)

        kx(ii,ind) = k0(ii)*sin(theta_i(ind)*pi/180);
        kz_air(ii,ind) = sqrt(k0(ii)^2-kx(ii,ind)^2);
        kz_end(ii,ind) = sqrt(k0(ii)^2*ep_allO(end)-kx(ii,ind)^2*ep_allO(end)/ep_allE(end));
        kz_EMT(ii,ind) = sqrt(k0(ii)^2*ep_iEMTO(ii)-kx(ii,ind)^2*ep_iEMTO(ii)/ep_iEMTE(ii));

%         for kk = 1:length(ep_allO)
%             kz_layer(ii,kk) = sqrt(k0(ii)^2*ep_allO(kk)-kx(ii,ind)^2*ep_allO(kk)/ep_allE(kk));
%             kz_norm(ii,kk) = kz_layer(ii,kk)/k0(ii);
%         end

        [A_TM{ii,ind}, B_TM{ii,ind}] = ABcoeff_aniso_multilayer(diff, layers, wl(ii), ep_allO, ep_allE, kx(ii,ind));
        [A_EMTp{ii,ind},B_EMTp{ii,ind}] = ABcoeff_aniso_multilayer(diff_EMT, 3, wl(ii), ep_3layerO, ep_3layerE, kx(ii,ind));
        for jj=1:length(z_norm)
            [~, ~, S_xTM(ii,jj), S_zTM(ii,jj)] = Poynting_aniso2(A_TM{ii}, B_TM{ii}, A_TM{ii}, B_TM{ii}, z_norm(jj), diff, diff_norm, kx(ii, ind), ep_allO, ep_allE, wl(ii));
        end
        for kk = 1:length(z_normEMT)
            [~,~,S_xEMTp(ii,kk), S_zEMTp(ii,kk)] = Poynting_aniso2(A_EMTp{ii}, B_EMTp{ii}, A_EMTp{ii}, B_EMTp{ii}, z_normEMT(kk), diff_EMT, [0, 1], kx(ii, ind), ep_3layerO, ep_3layerE, wl(ii));
        end
        x_normTM_temp(ii,:) = Streamline(S_xTM(ii,:),S_zTM(ii,:),z_norm); %znorm - z/d
        x_normTM(ii,:) = (x_normTM_temp(ii,:)-x_normTM_temp(ii,origin)); %x/d
    %     [x_normTE(ii,:)] = Streamline(S_xTE(ii,:),S_zTE(ii,:),z_norm);
        x_normEMT_temp(ii,:) = Streamline(S_xEMTp(ii,:),S_zEMTp(ii,:),z_normEMT);
        x_normEMT(ii,:) = (x_normEMT_temp(ii,:)-x_normEMT_temp(ii,origEMT));
    %     theta_s(ii) = -180/pi*atan(S_x(ii,origin)/S_z(ii,origin));
    %     theta_layer(ii) = -180/pi*atan((x_norm(ii,backside)-x_norm(ii,origin))/(z_norm(backside)-z_norm(origin)));
        
        Rf_TM(ii,ind) = abs(B_TM{ii,ind}(1)/A_TM{ii,ind}(1))^2;
        Tr_TM(ii,ind) = real(kz_end(ii,ind)/ep_allO(end))/real(kz_air(ii,ind)/ep_allO(1))*abs(A_TM{ii,ind}(end)/A_TM{ii,ind}(1))^2;
        Ab_TM(ii,ind) = 1-Rf_TM(ii,ind)-Tr_TM(ii,ind);

        Rf_EMT(ii,ind) = abs(B_EMTp{ii,ind}(1)/A_EMTp{ii,ind}(1))^2;
        Tr_EMT(ii,ind) = abs(A_EMTp{ii,ind}(end)/A_EMTp{ii,ind}(1))^2;
        Ab_EMT(ii,ind) = 1-Rf_EMT(ii,ind)-Tr_EMT(ii,ind);

%         depth_EMT(ii) = 1/(2*imag(kz_EMT(ii)))*1e6;
        
    end
    [Abs_EMT(ii),Abs_TM(ii)] = Hemispherical(Ab_EMT(ii,:),Ab_TM(ii,:),theta_i,d_theta);
    Abs_avg(ii) = (Abs_EMT(ii)+Abs_TM(ii))/2;
end
%[abs,ems] = Spectral(wl,Abs_TM); %need???




% figure (1); semilogx(wl,depth_EMT);
% figure (2); semilogx(wl,Tr_TM,'b',wl,Tr_EMT,'r--');
% figure (3); semilogx(wl,Rf_TM,'b',wl,Rf_EMT,'r--');
% figure (4); semilogx(wl,Ab_TM,'b',wl,Ab_EMT,'r--');
% figure (5); semilogx(wl,ep_iEMTO,'b',wl,ep_iEMTE,'r--');

%set(0,'defaulttextinterpreter','latex')
%figure (1)
%[Wave,Angle]=meshgrid(wl*1e3,theta_i);
% hand1=pcolor(Angle',Wave',log(R_123_s)/log(10)); %dont need
% set(hand1,'LineStyle','none');                   %dont need
%contourf(Angle',Wave',Rf_TM);
%xlabel('Incident Angle, $\theta_i$');
%ylabel('Wavelength, $\lambda$ (nm)');
%colormap hot;

jieshu=clock;
disp(datestr(jieshu));
calctime=etime(jieshu,kaishi);

toc;