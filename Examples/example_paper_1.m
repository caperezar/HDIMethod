% Convergence of the Harmonic Density Interpolation (HDI) method for 
% evaluation of on-surface Green's formulae

% Written by Carlos Perez-Arancibia (caperezar@gmail.com), Luiz Faria and
% Catalin Turc.

clc
close all
clear

%-------------------------------------------------------------------------%
% Definition of an exact Laplace solution and its gradient

% location of the source point lying outside the obstacle
xc = [2 2 2]; 

% computes the distance |x-xc|
R =  @(x,xc) vecnorm(x-repmat(xc,size(x,1),1),2,2);

% exact solution of the Laplace equation in the interior of the obstacle
F =  @(x,xc) 1./R(x,xc);

% components of the gradient of F
Fx = @(x,xc) -((x(:,1)-repmat(xc(1),size(x,1),1)))./R(x,xc).^3;
Fy = @(x,xc) -((x(:,2)-repmat(xc(2),size(x,1),1)))./R(x,xc).^3;
Fz = @(x,xc) -((x(:,3)-repmat(xc(3),size(x,1),1)))./R(x,xc).^3;


%-------------------------------------------------------------------------%
% Evaluation of Green's formula for various surface discretizations

% number of discretization points to be used
Narray = round(2.^(2:0.5:6.5)); 

% initialization of arrays that will contain the errors
ErrorSLDL = zeros(size(Narray)); 
ErrorDSHS = zeros(size(Narray));

for k=1:size(Narray,2)
    
    % each surface patch is discretized using a N x N Chebyshev grid
    N = Narray(k); 
    
    % contruction of the obstacle
    bdy = Body;
    bdy = bdy.cctor(N,'Cube',[],[]);
%     bdy = bdy.cctor(N,'Bean',[],[]);

    %---------------------------------------------------------------------%
    % contruction of the densities corresponding the exact Dirichlet (f) and
    % Neumann (fn) traces of the exact Laplace solution
    
    f = Density;
    f = f.cctor(bdy);

    fn = Density;
    fn = fn.cctor(bdy); 
   
        
    for j=1:bdy.nPat % loop over the surface patches
        
        % Dirichlet trace
        V = F(bdy.Xp(:,:,j),xc)-F(bdy.Xp(:,:,j),-xc);

        % Neumann trace
        Vx = Fx(bdy.Xp(:,:,j),xc)-Fx(bdy.Xp(:,:,j),-xc);
        Vy = Fy(bdy.Xp(:,:,j),xc)-Fy(bdy.Xp(:,:,j),-xc);
        Vz = Fz(bdy.Xp(:,:,j),xc)-Fz(bdy.Xp(:,:,j),-xc); 

        Vn = bdy.Nrm(:,1,j).*Vx + bdy.Nrm(:,2,j).*Vy + bdy.Nrm(:,3,j).*Vz;

        % fill the vector densities
        f  = f.matrix_to_density(j,V);   
        fn = fn.matrix_to_density(j,Vn);

    end
    
    % plot of the Diriclet trace 
    if k==6
        f.plot(bdy,1,'real')
    end
    
    %---------------------------------------------------------------------%
    % Evaluate all four boundary integral operators
    
    SLf = SL(fn,bdy); % single-layer 
    DLf = DL(f,bdy);  % double-layer
    DSf = DS(fn,bdy); % adjoint double-layer
    HSf = HS(f,bdy);  % hypersingular operator
    
    % error in evaluation of the Laplace solution on the boundary
    ErrorSLDL(k) = norm(SLf.vec - DLf.vec-f.vec/2,inf)/norm(f.vec,inf); 
    
    % error in the evalaution of the normal derivative of the Laplace
    % solution
    ErrorDSHS(k) = norm(DSf.vec - HSf.vec-fn.vec/2,inf)./norm(fn.vec,inf);
    
end


E3 = ErrorSLDL(end)*Narray(end)^3*Narray.^(-3); % third-order slope
E2 = ErrorDSHS(end)*Narray(end)^2*Narray.^(-2); % second-order slope

% plot of the errors in log-log scale
figure(2)
loglog(Narray,ErrorSLDL,'x-',Narray,ErrorDSHS,'x-',Narray,E3,'s-.',Narray,E2,'o-.')





