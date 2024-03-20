function [V,omegas,y,U,D1] = eig_OSS_temporal(ydata,Udata,alpha,beta,Re,Ncheb)
    % 3D wavenumber
    k2 = alpha^2+beta^2;

    % y bounds
    maxyd = max(ydata);
    minyd = min(ydata);  
    H = maxyd-minyd;

    % Derivatives
    [y,DM] = chebdif(Ncheb,3);
    % Derivatives without bcs
    D1 = DM(:,:,1); 
    D2 = DM(:,:,2);

    % 4th derivative with clamped bcs
    [~,D4c] = cheb4c(Ncheb); 

    % Scaling and interpolation
    scale = H/2;
    y = y*scale+(minyd-min(y*scale));
    U = interp1(ydata,Udata,y,'pchip',Udata(end));
    D1 = D1/scale;
    D2 = D2/(scale^2);
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
    D2c = D2(2:Ncheb-1,2:Ncheb-1);
    Ubc = U(2:Ncheb-1);
    Uy = Uy(2:Ncheb-1);
    Uyy = Uyy(2:Ncheb-1);

    II = eye(size(D2c));
    ZZ = zeros(size(D2c));

    Los = 1i*alpha*diag(Ubc)*(k2*II-D2c) ...
        + 1i*alpha*diag(Uyy) ...
        + 1/Re*(k2^2*II-2*k2*D2c+D4c);
    Lsq = 1i*alpha*diag(Ubc) ...
        + 1/Re*(k2*II-D2c);
   
    L = [Los,              ZZ; 
         1i*beta*diag(Uy), Lsq];
    M = [k2*II-D2c, ZZ; 
         ZZ,        II]; 

    [V,omegas]=eig(L,1i*M);
    omegas = diag(omegas);
    nV1 = size(V,1)/2;
    nV2 = size(V,2);
   
    % Pad ends with zeros
    V = [zeros(1,nV2); V(1:nV1,:); zeros(1,nV2); 
         zeros(1,nV2); V(nV1+1:end,:); zeros(1,nV2)];
end
