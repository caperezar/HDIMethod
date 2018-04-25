% Solution of a Laplace Neumann problem in the exterior of two touching 
% obstacles using the Harmonic Density Interpolation method

% Written by Carlos Perez-Arancibia (caperezar@gmail.com), Luiz Faria and
% Catalin Turc.

clear
clc
close all

% obstacle types:
type = ["Sphere","Cushion"];
% number of obstacles       
Nobs = numel(type);

% initialize array of obstacles
bdy(1:Nobs) = Body; 

% number of points used in the obstacle discretization (each patch is 
% discretized using N x N Chebyshev points)
N = 20;

% centers of the obstacles
center = [-0.103917010029620  -0.991950681940547  -0.077980529349595;
           0.396082989970380   1.462423318059453   0.511613034055861];
   
%-------------------------------------------------------------------------%
% Definition of an exact Laplace solution and its gradient

% point source locations
xc = center+[0.64 -0.18 0.03;0 0 0];

R =  @(x,y,z,o) sqrt((x-xc(o,1)).^2+(y-xc(o,2)).^2+(z-xc(o,3)).^2);
F =  @(x,y,z,o) 1./R(x,y,z,o);
Fx = @(x,y,z,o) -((x-xc(o,1)))./R(x,y,z,o).^3;
Fy = @(x,y,z,o) -((y-xc(o,2)))./R(x,y,z,o).^3;
Fz = @(x,y,z,o) -((z-xc(o,3)))./R(x,y,z,o).^3;

H = @(x,y,z) F(x,y,z,1)+F(x,y,z,2);
Hx = @(x,y,z) Fx(x,y,z,1)+Fx(x,y,z,2);
Hy = @(x,y,z) Fy(x,y,z,1)+Fy(x,y,z,2);
Hz = @(x,y,z) Fz(x,y,z,1)+Fz(x,y,z,2);


%-------------------------------------------------------------------------%
% initialize array of densities that will contain the exact Dirichlet and
% Neumann traces
f(1:Nobs)   = Density;
fn(1:Nobs)  = Density;

for o = 1:Nobs % loop over the obstacles
    
    % construct obstacle
    bdy(o)  = bdy(o).cctor(N,type(o),center(o,:));

    % construct density
    f(o)    = f(o).cctor(bdy(o));
    fn(o)   = fn(o).cctor(bdy(o));
        
    for j=1:bdy(o).nPat % loop over the patches

        % evaluate exact Dirichlet trace
        V = H(bdy(o).Xp(:,1,j),bdy(o).Xp(:,2,j),bdy(o).Xp(:,3,j));

        % evaluate exact Neumann trace
        Vx = Hx(bdy(o).Xp(:,1,j),bdy(o).Xp(:,2,j),bdy(o).Xp(:,3,j));
        Vy = Hy(bdy(o).Xp(:,1,j),bdy(o).Xp(:,2,j),bdy(o).Xp(:,3,j));
        Vz = Hz(bdy(o).Xp(:,1,j),bdy(o).Xp(:,2,j),bdy(o).Xp(:,3,j)); 

        Vn = bdy(o).Nrm(:,1,j).*Vx + bdy(o).Nrm(:,2,j).*Vy + bdy(o).Nrm(:,3,j).*Vz;

        % save into densities
        f(o)  = f(o).matrix_to_density(j,V);   
        fn(o) = fn(o).matrix_to_density(j,Vn);

    end

end


%-------------------------------------------------------------------------%
% Evaluation of the right-hand-side of the integral equation
SLfn = SL(fn,bdy);


%-------------------------------------------------------------------------%
% Iterative solution of the linear system
[dens,it] = GMRESsolve(SLfn,bdy);
% dens is an approximation of f

%-------------------------------------------------------------------------%
% Evaluate the approximate solution on planes

% construction of the planes
pl1 = Plane;
pl1 = pl1.cctor([0 0 0],[1 0 0;0 0 1],[4 4],[100 100]);

pl2 = Plane;
pl2=pl2.cctor([0 0 0],[0 1 0;0 0 1],[8 4],[200 100]);

% potentials computed with the HDI method
fld1  = evalDL(pl1.Xp,dens,bdy)-evalSL(pl1.Xp,fn,bdy);
fld2  = evalDL(pl2.Xp,dens,bdy)-evalSL(pl2.Xp,fn,bdy);

% potentials computed without using the HDI method
fld1p = evalDL2(pl1.Xp,dens,bdy)-evalSL2(pl1.Xp,fn,bdy);
fld2p = evalDL2(pl2.Xp,dens,bdy)-evalSL2(pl2.Xp,fn,bdy);

%-------------------------------------------------------------------------%
% Error evaluation

% errors on plane 1
error_fld1  = fld1 - H(pl1.Xp(:,1),pl1.Xp(:,2),pl1.Xp(:,3));
error_fld1p = fld1p - H(pl1.Xp(:,1),pl1.Xp(:,2),pl1.Xp(:,3));

% maximum errors
ee  = norm(error_fld1,inf)/norm(H(pl1.Xp(:,1),pl1.Xp(:,2),pl1.Xp(:,3)),inf);
eep  = norm(error_fld1p,inf)/norm(H(pl1.Xp(:,1),pl1.Xp(:,2),pl1.Xp(:,3)),inf);

% errors on plane 2
error_fld2  = fld2 - H(pl2.Xp(:,1),pl2.Xp(:,2),pl2.Xp(:,3));
error_fld2p = fld2p - H(pl2.Xp(:,1),pl2.Xp(:,2),pl2.Xp(:,3));


% surface errors on each obstacle
eDens(1) = norm(dens(1).vec-f(1).vec,inf)/norm(f(1).vec,inf);
eDens(2) = norm(dens(2).vec-f(2).vec,inf)/norm(f(2).vec,inf);
    

%-------------------------------------------------------------------------%
% Error plots

% plotting errors on planes
pl1.plot(error_fld1,1,'log10');
pl2.plot(error_fld2,1,'log10');
bdy(1).plot(1);bdy(2).plot(1);
caxis([-7 -0])
alpha  0.8;

pl1.plot(error_fld1p,2,'log10');
pl2.plot(error_fld2p,2,'log10');
bdy(1).plot(2);bdy(2).plot(2);
caxis([-7 -0])
alpha  0.8;



    