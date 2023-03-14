% clear all;
% close all;
addpath(genpath('Codes'));

%% Input parameters
H = 25;
Re = 1000;
beta = 0;
alpha = 0.3;
Ncheb = 256;

% Grid
[y,U,~,~] = blasius_profile(H,Ncheb);


%% Eigenproblem
[eigvects,eigvals,yinterp,Uinterp,D1] = eig_OSS_temporal(y,U,alpha,beta,Re,Ncheb);
sqidx = all(eigvects(1:length(y),:)==0);
sqvects = eigvects(:,sqidx);
sqvals = eigvals(sqidx);
osvects = eigvects(:,~sqidx);
osvals = eigvals(~sqidx);

figure()
h = plot(real(sqvals)/alpha,imag(sqvals)/alpha,'x', ...
         real(osvals)/alpha,imag(osvals)/alpha,'o'); 
grid on;
yline(0,'k-')
title(['Eigspectrum (Re = ', num2str(Re), ', \alpha = ',num2str(alpha), ', \beta = ',num2str(beta),', c = \omega/\alpha)']);
xlabel('$c_r$','Interpreter','latex'); ylabel('$c_i$','Interpreter','latex');
ylim([-2,2]);
xlim([-2,2]);
drawnow();
set(h,'ButtonDownFcn',{@plotTemporalEigenModesFromFigure,yinterp,Uinterp,D1,eigvects,eigvals,alpha,beta});