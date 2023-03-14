function plotSpatialEigenModesFromFigure(hObj,event,y,U,D1,eigvects,eigvals,beta)
%  FIND NEAREST (X,Y) COORDINATE TO MOUSE CLICK AND PLOT CORRESPONDING STABILITY MODE
% Inputs:
%  hObj (unused) the axes
%  event: info about mouse click
% OUTPUT
%  coordinateSelected: the (x,y,z) coordinate you selected
%  minIDx: The index of your inputs that match coordinateSelected. 

ny = length(y);

ar = hObj.XData; % Real part
ai = hObj.YData; % Imaginary part

pt = event.IntersectionPoint;       % The (x0,y0,z0) coordinate you just selected
pt = pt(1:end-1);                   % Remove z0
coordinates = [ar(:),ai(:)];        % matrix of your input coordinates
% distance between your selection and all points
dist = sqrt((pt(1)-coordinates(:,1)).^2+(pt(2)-coordinates(:,2)).^2);      
[~, minIdx] = min(dist);            % index of minimum distance to points

% Get index in eigenvalues list
alpha = coordinates(minIdx,1)+1i*coordinates(minIdx,2);
[~,eigIdx] = min(abs(eigvals-alpha));
vely = eigvects(1:ny,eigIdx);
eta = eigvects(ny+1:end,eigIdx);

% Dan's book page 42
k2 = alpha^2+beta^2;
velx = 1i/k2*(alpha*D1*vely-beta*eta);
velz = 1i/k2*(beta*D1*vely+alpha*eta);

% Reconstruction
xv = linspace(-2*pi/abs(real(alpha)),...
              2*pi/abs(real(alpha)),100);
vr = real(vely*exp(1i*alpha*xv));
ur = real(velx*exp(1i*alpha*xv));
wr = real(velz*exp(1i*alpha*xv));

%% FIGURE
f1 = figure(155);

% Plot eigenmode
subplot(1,5,1);
plot(real(vely),y,imag(vely),y,abs(vely),y);
legend('Re$(v)$','Im$(v)$','Abs$(v)$',"Interpreter","latex")
title(sprintf('$\\alpha= %.3f + %.3f i$',real(alpha),imag(alpha)),"Interpreter","latex");
ylabel('$y$',"Interpreter","latex")

subplot(1,5,2);
plot(real(eta),y,imag(eta),y,abs(eta),y);
legend('Re$(\eta)$','Im$(\eta)$','Abs$(\eta)$',"Interpreter","latex")
title(sprintf('$\\alpha= %.3f + %.3f i$',real(alpha),imag(alpha)),"Interpreter","latex");

% Plot u
subplot(1,5,3); hold on
pcolor(xv,y,ur); shading interp;
plot((U-(max(U)+min(U))/2)/(max(U)-min(U))*0.8*(max(xv)-min(xv))+(max(xv)+min(xv))/2,y,'k-','linewidth',2)
title('$u$',"Interpreter","latex")
xlabel('$x$',"Interpreter","latex");
xlim([min(xv),max(xv)]);
ylim([min(y),max(y)]);
if max(abs(ur(:))) ~= 0
    clim([-max(abs(ur(:))),max(abs(ur(:)))])
end
colormap('bluewhitered'); colorbar; 
hold off

% Plot v
subplot(1,5,4); hold on
pcolor(xv,y,vr); shading interp;
plot((U-(max(U)+min(U))/2)/(max(U)-min(U))*0.8*(max(xv)-min(xv))+(max(xv)+min(xv))/2,y,'k-','linewidth',2)
title('$v$',"Interpreter","latex")
xlabel('$x$',"Interpreter","latex");
xlim([min(xv),max(xv)]);
ylim([min(y),max(y)]);
if max(abs(vr(:))) ~= 0
    clim([-max(abs(vr(:))),max(abs(vr(:)))])
end
colormap('bluewhitered'); colorbar;
hold off

% Plot w
subplot(1,5,5); hold on
pcolor(xv,y,wr); shading interp;
plot((U-(max(U)+min(U))/2)/(max(U)-min(U))*0.8*(max(xv)-min(xv))+(max(xv)+min(xv))/2,y,'k-','linewidth',2)
title('$w$',"Interpreter","latex")
xlabel('$x$',"Interpreter","latex");
xlim([min(xv),max(xv)]);
ylim([min(y),max(y)]);
if max(abs(wr(:))) ~= 0
    clim([-max(abs(wr(:))),max(abs(wr(:)))])
end
colormap('bluewhitered'); colorbar; 
hold off

f1.Position = [ 900 80 1500 500];
end