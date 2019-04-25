% The script unpacks *.mat file (that was generated in FDTD) and makes
% contour plot of H and E fields.
% A.Yorzh 11/27/2018

% Load *.mat file
filename = "Bi2Se3_1400_diel1st.mat";
F = load(filename);

% Squeeze 3 dimensional arrays into 2-dimensional arrays
% squeeze - removes singleton dimensions
Ex = squeeze(F.Ex);
Ey = squeeze(F.Ey);
Ez = squeeze(F.Ez);

Hx = squeeze(F.Hx);
Hy = squeeze(F.Hy);
Hz = squeeze(F.Hz);

y = F.y - 100*1e-9;
z = F.z - 100*1e-9;

% Create an array of magnitudes of E and H
for i =1: size(Hx,1)*size(Hx,2)
    E(i) = norm([Ex(i), Ey(i), Ez(i)]);
    H(i) = norm([Hx(i), Hy(i), Hz(i)]);
end

% Now we have two big 1-dimensional arrays. We need to reshape them to the
% original shapes (check Ex,Ey,Ez, x, y, etc)
E = reshape(E, [size(Ex,1), size(Ex,2)]);
H = reshape(H, [size(Hx,1), size(Hx,2)]);

E = E./max(E(:));
H = H./max(H(:));

fig = figure;

% Here we collect data from Bi2Se3__MLTI_TRA_contour and put streamlines on
% fdtd plot.
[z_norm, x_normTM, depth] = Bi2Se3__MLTI_TRA_contour();
%[z_norm, x_normTM, depth] = Bi2Te3__MLTI_TRA_contour();
Pd = strcat('\delta =', num2str(depth,3),'microns');
angle_1 = ((90*pi)./180);
[theta, r] = cart2pol(x_normTM, z_norm);
theta = theta + angle_1;
[x_norm, z_norm] = pol2cart(theta, r);
x_norm = x_norm*2.0*1E-6;
z_norm = z_norm*2.0*1E-6;


% Make contour plot for E
hold off
subplot(2,1,1)
c1 = contourf(z, y, E,'edgecolor','none');
%title('E Bi2Se3 0.5um')
set(gca, 'Ydir', 'reverse');
ylim([-1E-5 1E-5]);
xticks([-2E-6 -1E-6 0 1E-6 2E-6 3E-6])
%text(2.5E-6, -9.5E-6, '(a)', 'Color', 'white', 'Fontsize', 16)
xticklabels({'-2', '-1.0','0.0', '1.0', '2.0', '3.0'})
yticklabels({'-10', '-5', '0.0', '5', '10'})
xlabel ('z (microns)')
%ylabel ('y (um)')
view(-90,90);
colormap('jet');
%set(colorbar, 'TickLabels', {'0', '1','2', '3', '4', '5', '6'} );
%colorbar;


% Draw a white line to show the slab
hold on
a1 = [0.0 0.0];
b1 = [-1.5e-5 1.5e-5];
plot (a1,b1, 'color','w', 'linewidth', 1.2)
a2 = [-1.98e-6 -1.98e-6];
b2 = [-1.5e-5 1.5e-5];
plot (a2,b2, 'color','w', 'linewidth', 1.5)
pbaspect([1 2.5 1])
text(2.5E-6, 7.7E-6, '||E||', 'Color', 'white', 'Fontsize', 16)
%text(-1.5E-6, 3.7E-6, Pd, 'Color', 'white', 'Fontsize', 12)
% Plot streamline on contour plot
plot(x_norm, z_norm, 'Color', 'black', 'linewidth', 2, 'linestyle', '--');
xlim([-2E-6, 3E-6]);

% Make contour plot for H
hold off
subplot(2,1,2)
c2 = contourf(z, y, H,'edgecolor','none');
ylim([-1.0E-5 1.0E-5]);
xticks([-2E-6 -1E-6 0 1E-6 2E-6 3E-6])
xticklabels({'-2', '-1.0','0.0', '1.0', '2.0', '3.0'})
yticklabels({'-10', '-5', '0.0', '5', '10'})
xlabel ('z (microns)')
ylabel ('y (microns)')
%title('H Bi2Se3 0.5um')
%title('H Bi2Te3 0.962um')
set(gca, 'Ydir', 'reverse');
view(-90,90);
colormap('jet')
%colorbar;
pbaspect([1 2.5 1])

% Draw a white line to show the slab
hold on
a1 = [0.0 0.0];
b1 = [-1.5e-5 1.5e-5];
%set(gca, 'XTickLabel', get(gca, 'XTick'));
plot (a1,b1, 'color','w', 'linewidth', 1.2)
a2 = [-1.98e-6 -1.98e-6];
b2 = [-1.5e-5 1.5e-5];
plot (a2,b2, 'color','w', 'linewidth', 1.5)
%text(-1.5E-6, 3.7E-6, Pd, 'Color', 'white', 'Fontsize', 12)
%text(2.5E-6, -9.5E-6, '(b)', 'Color', 'white', 'Fontsize', 16)
text(2.5E-6, 7.7E-6, '||H||', 'Color', 'white', 'Fontsize', 16)
% Plot streamline on contour plot
plot(x_norm, z_norm, 'Color', 'black', 'linewidth', 2, 'linestyle', '--');
xlim([-2E-6, 3E-6]);
hold off



%test
%figure;
%plot(x_normTM*2.0*1E-6, z_norm*2.0*1E-6);
%set(gca, 'Ydir', 'reverse');
