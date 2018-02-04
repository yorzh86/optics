% Author: Richard Zihao Zhang (03/02/15)
% Calculation of Fabry-Perot interferometer using both surface current
% reflectance, or multilayer SMM.
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

% wl = 0.6, 10 and 20 um (first two from Lehman Nano Lett)
wl = [0.4];
nu = 1e4./wl;
%om = 2*pi*c0./wl.*1e6;
%k0 = om./c0;
theta_i = [0, 10, 20, 30];       %degrees

per = 20;
graphene = 1;
d_g = graphene*0.335e-9;
d_d = 100e-9;
% d_g = 20e-9;
% d_d = 20e-9;
diff = zeros(1,2*per+2);
diff(2:2:end) = d_g;
diff(3:2:end-1) = d_d;

% bigstep = 500;
% smallstep = 50;
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
dif_n = [0, 1];

origEMT = find(z_normEMT==0,1,'first');

ff = d_g/(d_g+d_d);

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
T1 = 300;
% mu=0.3633;  % enz @ 15 um
mu=0.4157;

%-DIELECTRIC EPSILON-----------------------------------
ep_SiO2 = 2.2;
ep_dO = ep_SiO2;
ep_dE = ep_SiO2;

% [ep_hBNO, ep_hBNE] = epsilon_hBN(om);
% ep_dO = ep_hBNO;
% ep_dE = ep_hBNE;

%-METAL EPSILON-----------------------------------
%[sigma,~]=Graphene_conductivity(om,T1,mu);
% ep_g=1i*sigma/(om*ep0*d_g);
% % ep_Ag = epsilon_Ag(om);
% ep_mO = ep_g+ep_dO;
% ep_mE = ep_dE;

% if (wl>9.9)
%   ep_Au = epsi_Au_Drude(om);
% else
%   % Au up to 10 um
%   ep_Au = Palik_Au(om);
% end
% ep_Si = epsi_Si(om, 300, 0, 0);

% ep_pO = [ep_mO, ep_dO];
% ep_pE = [ep_mE, ep_dE];
ep1 = 1.0;
ep3 = 1.0;

%------------------------------------------

% ep_MLO = repmat(ep_pO,1,per);
% ep_MLE = repmat(ep_pE,1,per);
% ep_O = [ep1,ep_MLO,ep_mO,ep3];
% ep_E = [ep1,ep_MLE,ep_mE,ep3];
% layers = length(ep_O);
% 
% ep_EMTO = ff*ep_mO+(1-ff)*ep_dO;
% ep_EMTE = ep_mE*ep_dE/(ff*ep_dE+(1-ff)*ep_mE);
% 
% ep_3O = [ep1, ep_EMTO, ep3];
% ep_3E = [ep1, ep_EMTE, ep3];
%------------------------------------------
diff = [1E6, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03,...
    0.03, 0.03, 0.03, 0.03, 0.03, 1E6];
layers = 18;
%ep_mO = [3, 3, 3, 4, 1];
%ep_mE = [2, 2, 1, 1, 1];
ep_mO =[-7.7123+0.0505i, -12.6127+0.0986i,-18.6019+0.1705i, -25.6796+0.2707i,-33.8458+0.4040i];
ep_mE = [-7.7123+0.0505i, -12.6127+0.0986i,-18.6019+0.1705i, -25.6796+0.2707i,-33.8458+0.4040i];

ep_dO =[9.1038+0.0782i, 7.3514,  6.7857, 6.5088, 6.3491];
ep_dE =[11.4552+0.5068i, 9.2018, 8.4019, 8.0158, 7.7948];
%ep_dO = [1, 2, 3, 4, 5];
%ep_dE = [2, 2, 3, 3, 4];

for ii = 1:length(theta_i)
    for jj = 1:length(wl)
        ep_O = [ep1, ep_mO(jj), ep_dO(jj), ep_mO(jj),ep_dO(jj), ep_mO(jj),ep_dO(jj), ep_mO(jj),...
            ep_dO(jj), ep_mO(jj),ep_dO(jj), ep_mO(jj),ep_dO(jj), ep_mO(jj),...
            ep_dO(jj), ep_mO(jj),ep_dO(jj), ep3];
        
        ep_E = [ep1, ep_mE(jj), ep_dE(jj),ep_mE(jj), ep_dE(jj),ep_mE(jj), ep_dE(jj),ep_mE(jj), ep_dE(jj),...
            ep_mE(jj), ep_dE(jj),ep_mE(jj), ep_dE(jj),ep_mE(jj), ep_dE(jj),ep_mE(jj), ep_dE(jj),...
            ep3];
        om = 2*pi*c0./wl(jj).*1e6;
        k0 = om./c0;
    
        kx(ii) = k0*sin(theta_i(ii)*pi/180);
        kz_air(ii) = sqrt(k0^2-kx(ii)^2);
        kz_end(ii) = sqrt(k0^2*ep_O(end)-kx(ii)^2*ep_O(end)/ep_E(end));
        %for kk = 1:length(ep_O)
        %    kz_layer(ii,kk) = sqrt(k0^2*ep_O(kk)-kx(ii)^2*ep_O(kk)/ep_E(kk));
        %    kz_norm(ii,kk) = kz_layer(ii,kk)/k0;
        %end

        [A_TM{ii}, B_TM{ii}] = ABcoeff_aniso_multilayer(diff, layers, wl(jj), ep_O, ep_E, kx(ii));
    
    %[A_EMTp{ii},B_EMTp{ii}] = ABcoeff_aniso_multilayer(diff_EMT, 3, wl,
    %ep_3O, ep_3E, kx(ii));1

%     for jj=1:length(z_norm)
%         [~, ~, S_xTM(ii,jj), S_zTM(ii,jj)] = Poynting_aniso2(A_TM{ii}, B_TM{ii}, A_TM{ii}, B_TM{ii}, z_norm(jj), diff, diff_norm, kx(ii), ep_O, ep_E, wl);
%     end
%     for kk = 1:length(z_normEMT)
%         [~,~,S_xEMTp(ii,kk), S_zEMTp(ii,kk)] = Poynting_aniso2(A_EMTp{ii}, B_EMTp{ii}, A_EMTp{ii}, B_EMTp{ii}, z_normEMT(kk), diff_EMT, [0, 1], kx(ii), ep_3O, ep_3E, wl);
%     end
    
%     x_normTM_temp(ii,:) = Streamline(S_xTM(ii,:),S_zTM(ii,:),z_norm);
%     x_normTM(ii,:) = (x_normTM_temp(ii,:)-x_normTM_temp(ii,origin));
% %     [x_normTE(ii,:)] = Streamline(S_xTE(ii,:),S_zTE(ii,:),z_norm);
%     x_normEMT_temp(ii,:) = Streamline(S_xEMTp(ii,:),S_zEMTp(ii,:),z_normEMT);
%     x_normEMT(ii,:) = (x_normEMT_temp(ii,:)-x_normEMT_temp(ii,origEMT));
    
%     theta_s(ii) = -180/pi*atan(S_x(ii,origin)/S_z(ii,origin));
%     theta_layer(ii) = -180/pi*atan((x_norm(ii,backside)-x_norm(ii,origin))/(z_norm(backside)-z_norm(origin)));
        Rf_TM(ii) = abs(B_TM{ii}(1)/A_TM{ii}(1))^2;
        Tr_TM(ii) = real(kz_end(ii)/ep_O(end))/real(kz_air(ii)/ep_O(1))*abs(A_TM{ii}(end)/A_TM{ii}(1))^2;
        Ab_TM(ii) = 1-Rf_TM(ii)-Tr_TM(ii);
    end
    
%     Rf_EMT(ii) = abs(B_EMTp{ii}(1)/A_EMTp{ii}(1))^2;
%     Tr_EMT(ii) = abs(A_EMTp{ii}(end)/A_EMTp{ii}(1))^2;
%     Ab_EMT(ii) = 1-Rf_EMT(ii)-Tr_EMT(ii);

end

disp(Rf_TM)
disp(Tr_TM)
%disp(Ab_TM)

% xx_EMT = x_normEMT;
% xx = x_normTM;
% zz = z_norm;
% zz_EMT = z_normEMT;
% 
% 
% figure (2);
% plot(xx,zz,xx_EMT,zz_EMT,'o');
% xlabel('z/slab');
% ylabel('x/slab');


% col_grad = Tr_TM;
% figure (3);
% surface(z_norm,x_normTM,Tr_TM,col_grad,'facecol','no','edgecol','interp');
% xlabel('z/slab');
% ylabel('x/slab');
% colormap hot;
% colormap(flipud(colormap))

% jieshu=clock;
% disp(datestr(jieshu));
% calctime=etime(jieshu,kaishi);
% 
% toc;