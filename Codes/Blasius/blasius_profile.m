function [y,U,D,W] = blasius_profile(H,ny)
    % Blasius solutions
    [~,Fp,eta] = blasius_sol();

    % Generate velocity fields
    U0 = 1;
    K = trapz(eta,1-Fp); %Displacement thickness factor

    [D,y] = cheb(ny-1);
    y = (y+1)*H/2;
    D = D/(H/2);

    [~,W] = clencurt(ny-1);
    W = W*(H/2);

    y = flip(y);
    eta0 = K*y;
    fp = interpFp(eta,Fp,eta0);
    U = U0*flip(fp);
    y = flip(y);
end

function [D,x] = cheb(N)
    if N==0, D=0; x=1; return, end
    x = cos(pi*(0:N)/N)';
    c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
    X = repmat(x,1,N+1);
    dX = X-X';
    D = (c*(1./c)')./(dX+(eye(N+1))); % off-diagonal entries
    D = D - diag(sum(D')); % diagonal entries
end

function [x,w] = clencurt(N)
    theta = pi*(0:N)'/N; x = cos(theta);
    w = zeros(1,N+1); ii = 2:N; v = ones(N-1,1);
    if mod(N,2)==0
       w(1) = 1/(N^2-1); w(N+1) = w(1);
       for k=1:N/2-1, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
       v = v - cos(N*theta(ii))/(N^2-1);
    else
       w(1) = 1/N^2; w(N+1) = w(1);
       for k=1:(N-1)/2, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
    end
    w(ii) = 2*v/N;
    w = w';
end

function bl=blasius_eq(eta,F)
    % Blasius equation
    % f*fpp + fppp = 0
    % F = [f,fp,fpp]
    bl = [F(2);F(3);-0.5*F(1)*F(3)]; % [fp,fpp,fppp] 
end

function err=findfpp(x)
    % Function to find the initial condition of fpp 
    % that yields f(eta_max) = 1
    [~,yout]=ode45(@blasius_eq,[0 10],[0 0 x]);
    err = (yout(end,2)-1)^2;
end

function [F,Fp,eta] = blasius_sol()
    [eta,~]=chebdif(2048,1);
    eta = flip((eta+1)*25); 
    options = optimset('TolX',1e-10);
    fpp0 = fminsearch(@findfpp,0.332,options);
    F0 = [0;0;fpp0];
    [~,F]=ode45(@blasius_eq,eta,F0);
    Fp = F(:,2);
    F = F(:,1);
end

function val = interpFp(eta,Fp,x)
   [~,index] = min(abs(x-max(eta)));
   x1 = x(1:index-1);
   x2 = x(index:end);
   val1 = interp1(eta,Fp,x1);
   val2 = ones([length(x2),1])*Fp(end);
   val = [val1;val2];
end