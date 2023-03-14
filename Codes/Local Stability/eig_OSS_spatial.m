function [V,alphas,y,U,D1] = eig_OSS_spatial(ydata,Udata,omega,beta,Re,Ncheb)
    % y bounds
    maxyd = max(ydata);
    minyd = min(ydata);  
    H = maxyd-minyd;

    % Derivatives
    [y,DM] = chebdif(Ncheb,3);
    % Derivatives without bcs
    D1 = DM(:,:,1); 
    D2 = DM(:,:,2);
    D3 = DM(:,:,3);

    % 4th derivative with clamped bcs
    [~,D4c] = cheb4c(Ncheb); 

    % Scaling and interpolation
    scale = H/2;
    y = y*scale+(minyd-min(y*scale));
    U = interp1(ydata,Udata,y,'spline',Udata(end)); % C2 interpolation
    D1 = D1/scale;
    D2 = D2/(scale^2);
    D3 = D3/(scale^3);
    D4c = D4c/(scale^4);

    % Computing derivatives of base flow without bcs
    Uy = D1*U;
    Uyy = D2*U; 

    % DEBUG
%     figure(100); clf;
%     subplot(1,3,1)
%     plot(U,y,'o-'); title('U');
%     subplot(1,3,2)
%     plot(Uy,y,'o-'); title('Uy');
%     subplot(1,3,3)
%     plot(Uyy,y,'o-'); title('Uyy');

    % Cut-off boundaries
    D1c = D1(2:Ncheb-1,2:Ncheb-1);
    D2c = D2(2:Ncheb-1,2:Ncheb-1);
    D3c = D3(2:Ncheb-1,2:Ncheb-1);
    Ubc = U(2:Ncheb-1);
    ybc = y(2:Ncheb-1);
    Uy = Uy(2:Ncheb-1);
    Uyy = Uyy(2:Ncheb-1);

    [s1,s2] = size(D1c);
    I = eye([s1,s2]);
    Z = zeros([s1,s2]);
    
    % Based on Schmid and Henningson - 2001
    R0 = 1i*omega*(D2c-beta^2*I)+1/Re*(D4c-2*beta^2*D2c+beta^4*I);
    R1 = -2i*omega*D1c - 4/Re*(D3c-beta^2*D1c) - 1i*diag(Ubc)*(D2c-beta^2*I) + 1i*diag(Uyy); 
    R2 = 4/Re*D2c + 2i*diag(Ubc)*D1c;

    T0 = -1i*omega*I - 1/Re*(D2c-beta^2*I);
    T1 = 2/Re*D1c + 1i*diag(Ubc);
    S = 1i*beta*diag(Uy);
    L = [-R1 -R0 Z; I Z Z; Z -S -T0];
    F = [R2 Z Z; Z I Z; Z Z T1];

    [V,alphas]=eig(L,F);
    alphas = diag(alphas);
    nV1 = size(V,1)/3;
    nV2 = size(V,2);
   
    %Haj-Hariri 1988 transformation
    V = [zeros(1,nV2); V(nV1+1:2*nV1,:).*exp(-ybc*alphas.'); zeros(1,nV2); 
         zeros(1,nV2); V(2*nV1+1:end,:).*exp(-ybc*alphas.'); zeros(1,nV2)]; % Pad ends with zeros
end