clear all;
% close all;
addpath(genpath('Codes'));

%% Input parameters
H = 20;
Re = 1000;
beta = 0;
omega = 0.26;
Ncheb = 256;

% Grid
[y,U,~,~] = blasius_profile(H,Ncheb);

%% Eigenproblem
[eigvects,eigvals,yinterp,Uinterp,D1] = eig_OSS_spatial(y,U,omega,beta,Re,Ncheb);
sqidx = all(eigvects(1:length(y),:)==0);
sqvects = eigvects(:,sqidx);
sqvals = eigvals(sqidx);
osvects = eigvects(:,~sqidx);
osvals = eigvals(~sqidx);

figure()
h = plot(real(sqvals),imag(sqvals),'x', ...
         real(osvals),imag(osvals),'o'); 
grid on;
yline(0,'k-')
title(['Eigspectrum (Re = ', num2str(Re), ', \beta = ',num2str(beta), ', \omega = ',num2str(omega),')']);
xlabel('\alpha_r'); ylabel('\alpha_i');
ylim([-1,1]);
xlim([-1,1]);
drawnow();
set(h,'ButtonDownFcn',{@plotSpatialEigenModesFromFigure,yinterp,Uinterp,D1,eigvects,eigvals,beta});

