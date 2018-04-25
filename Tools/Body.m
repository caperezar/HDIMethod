classdef Body
    
properties
    % number of points per grid per patch
    nPts
    
    % geometric variables at grid points
    Xp
    dXdu
    dXdv
    Nrm    
         
    %%
    h % parameter used in the finite difference discretizations
    
    % jacobian and integration weights
    Wp
    Jp
    
    % number of patches
    nPat
    
    % vectors that generate the Chebyshev grids
    tGrid
    chebGrid
    
    % paramtrization functions
    X
    dX
    d2X
    invX
    dn
    
    % total number of discretization points
    nPtot
    
    tol    
    %% Save geometry as densities
    save_opt
    
    dPu
    dPv
    
    d2Puu
    d2Puv
    d2Pvv
    
    d3Puuu
    d3Puuv
    d3Puvv
    d3Pvvv
    
    NR
    dNRu
    dNRv
        
    d2NRuu
    d2NRuv
    d2NRvv   
    
    d3NRuuu
    d3NRuuv
    d3NRuvv
    d3NRvvv
    
end

methods
    
    function obj=cctor(obj,n,name,center,rad,save_opt)
        
        srf= eval(name);
        
        obj.h = srf.h;
        
        obj.tol = 4*(2/n);
                
        if exist('center','var') && ~isempty(center)
                 
            srf.center = center;
            
        end
        
        if exist('rad','var') && ~isempty(rad) 
            srf.rad= rad;
        end
        
        obj.X  = @(u,v,id) srf.X(u,v,id);
        obj.dX = @(u,v,id) srf.dX(u,v,id);
        obj.d2X= @(u,v,id) srf.d2X(u,v,id);
        obj.invX= @(x,id) srf.invX(x,id);
        obj.dn= @(u,v,id) srf.dn(u,v,id);
        
        if size(n,1)==1
            
            obj.nPts = n*ones(6,1);
            
        else
            
            obj.nPts = n;
            
        end
        
        obj.nPat = srf.nPat;  
        
        [u,w,t] = fejer1(n);
        
        v=u;
        
        [U,V] = meshgrid(u,v);
        
        
        W = w*w';

        ut=reshape(U,n^2,1);
        vt=reshape(V,n^2,1);
        obj.Wp=reshape(W,n^2,1); % Weights of quadrature
        obj.tGrid = t;           % Uniform point in theta
        obj.chebGrid = u;        % Uniform point in theta
        
        
        for j=1:obj.nPat            
            obj.Xp(:,:,j) = srf.X(ut,vt,j); 
            [obj.dXdu(:,:,j),obj.dXdv(:,:,j),obj.Nrm(:,:,j),obj.Jp(:,j)] = srf.dX(ut,vt,j);            
        end
            
        if exist('save_opt','var') && save_opt==1
            
            obj.save_opt=1;
            
            P = Density;P = P.cctor(obj,3);
            P.vec = [];
            for j=1:obj.nPat            
                P.vec = [P.vec; srf.X(ut,vt,j)];                
            end
            
            obj.dPu = P.diff(obj,'U');
            obj.dPv = P.diff(obj,'V');

            e1_e2 = cross(obj.dPu.vec,obj.dPv.vec);   % Normal vector    
            jac = vecnorm(e1_e2,2,2) ; % Unit normal

            obj.NR= Density; obj.NR= obj.NR.cctor(obj,3);
            obj.NR.vec= e1_e2 ./ jac ; % Unit normal
            
            obj.d2Puu = obj.dPu.diff(obj,'U');
            obj.d2Puv = obj.dPu.diff(obj,'V');
            obj.d2Pvv = obj.dPv.diff(obj,'V');
            
            obj.d3Puuu = obj.d2Puu.diff(obj,'U');
            obj.d3Puuv = obj.d2Puv.diff(obj,'U');
            obj.d3Puvv = obj.d2Pvv.diff(obj,'U');
            obj.d3Pvvv = obj.d2Pvv.diff(obj,'V');

            obj.dNRu = obj.NR.diff(obj,'U');
            obj.dNRv = obj.NR.diff(obj,'V');
            
            obj.d2NRuu= obj.dNRu.diff(obj,'U');
            obj.d2NRuv= obj.dNRv.diff(obj,'U');
            obj.d2NRvv= obj.dNRv.diff(obj,'V');
            
            obj.d3NRuuu= obj.d2NRuu.diff(obj,'U');
            obj.d3NRuuv= obj.d2NRuv.diff(obj,'U');
            obj.d3NRuvv= obj.d2NRvv.diff(obj,'U');
            obj.d3NRvvv= obj.d2NRvv.diff(obj,'V');
            
        end
                
        
        %%%               
        obj.nPtot = sum(obj.nPts.^2);
    end
    
    
    function [n_uu,n_uv,n_vv] = d2n(obj,u,v,id)
        % second order derivatives of the normal
        h0 = obj.h;
%         h = 1.0e-5;

        [n2_u,n2_v] = obj.dn(u+h0,v,id);
        [n1_u,n1_v] = obj.dn(u-h0,v,id);
        
        n_uu = 0.5/h0*(n2_u-n1_u);
        n_uv = 0.5/h0*(n2_v-n1_v);

        [n2_u,n2_v] = obj.dn(u,v+h0,id);
        [n1_u,n1_v] = obj.dn(u,v-h0,id);

        n_uv = 0.5*n_uv + 0.5^2/h0*(n2_u-n1_u);
        n_vv = 0.5/h0*(n2_v-n1_v);

    end
    
    function [n_uuu,n_uuv,n_uvv,n_vvv]=d3n(obj,u,v,id)
        % third order derivatives of the normal
        h0 = obj.h;

        [n2_uu,n2_uv,n2_vv] = obj.d2n(u+h0,v,id);
        [n1_uu,n1_uv,n1_vv] = obj.d2n(u-h0,v,id);

        n_uuu = 0.5/h0*(n2_uu-n1_uu);
        n_uuv = 0.5/h0*(n2_uv-n1_uv);
        n_uvv = 0.5/h0*(n2_vv-n1_vv);

        [~,n2_uv,n2_vv]= obj.d2n(u,v+h0,id);
        [~,n1_uv,n1_vv]= obj.d2n(u,v-h0,id);

        n_uvv = 0.5*n_uvv + 0.5^2/h0*(n2_uv-n1_uv);
        n_vvv = 0.5/h0*(n2_vv-n1_vv);

    end 
    
    function [Xuuu,Xuuv,Xuvv,Xvvv] = d3X(obj,u,v,id)
       % third order derivatives of the parametrization
%        h = 1.0d-5;
       h0 = obj.h;

       [X2_uu,X2_uv,X2_vv] = obj.d2X(u+h0,v,id);
       [X1_uu,X1_uv,X1_vv] = obj.d2X(u-h0,v,id);

       Xuuu = 0.5/h0*(X2_uu-X1_uu);
       Xuuv = 0.5/h0*(X2_uv-X1_uv);
       Xuvv = 0.5/h0*(X2_vv-X1_vv);

      [~,X2_uv,X2_vv] = obj.d2X(u,v+h0,id);
      [~,X1_uv,X1_vv] = obj.d2X(u,v-h0,id);

      Xuvv = 0.5*Xuvv + 0.5^2/h0*(X2_uv-X1_uv);
      Xvvv = 0.5/h0*(X2_vv-X1_vv);

    end
    
    function [out,cA] = invAOld(obj,u,v,id) % HOSS interpolation matrix
        [dxdu,dxdv,n] = obj.dX(u,v,id);
        [dndu, dndv] = obj.dn(u,v,id);   
        [d2xdu2, d2xdudv, d2xdv2] = obj.d2X(u,v,id);
            
         A = [1, 0, 0, 0, 0, 0, 0, 0, 0;...
            0, dxdu(1), dxdu(2), dxdu(3), 0, 0, 0, 0, 0;...
            0, dxdv(1), dxdv(2), dxdv(3), 0, 0, 0, 0, 0;...
            0, n(1),n(2),n(3), 0, 0, 0, 0, 0;...
            ...
            0, dndu(1), dndu(2), dndu(3), ...
            n(1)*dxdu(2) +  n(2)*dxdu(1),...
            n(1)*dxdu(3) +  n(3)*dxdu(1),...
            n(2)*dxdu(3) +  n(3)*dxdu(2),...
            2*n(1)*dxdu(1) + 2*n(2)*dxdu(2) - 4*n(3)*dxdu(3),...
            2*n(1)*dxdu(1) - 2*n(2)*dxdu(2);...
            ...
            0, dndv(1), dndv(2), dndv(3), ...
            n(1)*dxdv(2) +  n(2)*dxdv(1),...
            n(1)*dxdv(3) +  n(3)*dxdv(1),...
            n(2)*dxdv(3) +  n(3)*dxdv(2),...
            2*n(1)*dxdv(1) + 2*n(2)*dxdv(2) - 4*n(3)*dxdv(3),...
            2*n(1)*dxdv(1) - 2*n(2)*dxdv(2);...
            ...
            0, d2xdu2(1),d2xdu2(2),d2xdu2(3),...
            2*dxdu(1)*dxdu(2),2*dxdu(1)*dxdu(3),2*dxdu(2)*dxdu(3),...
            2*(dxdu(1)^2 + dxdu(2)^2 - 2*dxdu(3)^2), 2*(dxdu(1)^2 - dxdu(2)^2);
            ...
            0, d2xdv2(1),d2xdv2(2),d2xdv2(3),...
            2*dxdv(1)*dxdv(2),2*dxdv(1)*dxdv(3),2*dxdv(2)*dxdv(3),...
            2*(dxdv(1)^2 + dxdv(2)^2 - 2*dxdv(3)^2), 2*(dxdv(1)^2 - dxdv(2)^2);...
            ...
            0, d2xdudv(1),d2xdudv(2),d2xdudv(3),...
            dxdu(1)*dxdv(2) + dxdv(1)*dxdu(2),...
            dxdu(1)*dxdv(3) + dxdv(1)*dxdu(3),...
            dxdu(2)*dxdv(3) + dxdv(2)*dxdu(3),...
            2*(dxdu(1)*dxdv(1) + dxdu(2)*dxdv(2) -2*dxdu(3)*dxdv(3)),...
            2*(dxdu(1)*dxdv(1) - dxdu(2)*dxdv(2))];
        
            cA = cond(A);            
            out = inv(A);                 
     end
%     
     function [out,cA,A] = invA(obj,p,q,id) % HDI interpolation matrix
         u = obj.chebGrid(p);
         v = obj.chebGrid(q);
         
        [dxdu,dxdv,n] = obj.dX(u,v,id);
        [dndu, dndv] = obj.dn(u,v,id);   
        [d2xdu2, d2xdudv, d2xdv2] = obj.d2X(u,v,id);
            
         A = [1, 0, 0, 0, 0, 0, 0, 0, 0;...
            0, dxdu(1), dxdu(2), dxdu(3), 0, 0, 0, 0, 0;...
            0, dxdv(1), dxdv(2), dxdv(3), 0, 0, 0, 0, 0;...
            0, n(1),n(2),n(3), 0, 0, 0, 0, 0;...
            ...
            0, dndu(1), dndu(2), dndu(3), ...
            n(1)*dxdu(2) +  n(2)*dxdu(1),...
            n(1)*dxdu(3) +  n(3)*dxdu(1),...
            n(2)*dxdu(3) +  n(3)*dxdu(2),...
            2*n(1)*dxdu(1) + 2*n(2)*dxdu(2) - 4*n(3)*dxdu(3),...
            2*n(1)*dxdu(1) - 2*n(2)*dxdu(2);...
            ...
            0, dndv(1), dndv(2), dndv(3), ...
            n(1)*dxdv(2) +  n(2)*dxdv(1),...
            n(1)*dxdv(3) +  n(3)*dxdv(1),...
            n(2)*dxdv(3) +  n(3)*dxdv(2),...
            2*n(1)*dxdv(1) + 2*n(2)*dxdv(2) - 4*n(3)*dxdv(3),...
            2*n(1)*dxdv(1) - 2*n(2)*dxdv(2);...
            ...
            0, d2xdu2(1),d2xdu2(2),d2xdu2(3),...
            2*dxdu(1)*dxdu(2),2*dxdu(1)*dxdu(3),2*dxdu(2)*dxdu(3),...
            2*(dxdu(1)^2 + dxdu(2)^2 - 2*dxdu(3)^2), 2*(dxdu(1)^2 - dxdu(2)^2);
            ...
            0, d2xdv2(1),d2xdv2(2),d2xdv2(3),...
            2*dxdv(1)*dxdv(2),2*dxdv(1)*dxdv(3),2*dxdv(2)*dxdv(3),...
            2*(dxdv(1)^2 + dxdv(2)^2 - 2*dxdv(3)^2), 2*(dxdv(1)^2 - dxdv(2)^2);...
            ...
            0, d2xdudv(1),d2xdudv(2),d2xdudv(3),...
            dxdu(1)*dxdv(2) + dxdv(1)*dxdu(2),...
            dxdu(1)*dxdv(3) + dxdv(1)*dxdu(3),...
            dxdu(2)*dxdv(3) + dxdv(2)*dxdu(3),...
            2*(dxdu(1)*dxdv(1) + dxdu(2)*dxdv(2) -2*dxdu(3)*dxdv(3)),...
            2*(dxdu(1)*dxdv(1) - dxdu(2)*dxdv(2))];
        
            cA = cond(A);            
            out = inv(A); 
     end
    
    
    function l = loc_to_glob(obj,jU,jV,id) 
        
        % i --> v
        % j --> u
        N0=0;
        
        if id>0
            
            N0= sum(obj.nPts(1:id-1).^2);
            
        end
        
        l = N0+jV+(jU-1)*obj.nPts(id);
    end
    
    
    function []=plot(obj,nfig)
        
        figure(nfig);hold on;
        
        for j=1:obj.nPat
            n = obj.nPts(j);
            X0=reshape(obj.Xp(:,1,j),n,n);
            Y0=reshape(obj.Xp(:,2,j),n,n);
            Z0=reshape(obj.Xp(:,3,j),n,n);
            mesh(X0,Y0,Z0,zeros(size(Z0)),'edgecolor','black');

        end
        axis tight
        axis equal
        hold off
    end
    
    function [I,II,Xu,Xv,Nrm]=geo_variables(obj,u,v,p)

        [Xu,Xv,Nrm] = obj.dX(u,v,p);
        [Xuu,Xuv,Xvv] = obj.d2X(u,v,p);

        E = dot(Xu,Xu,2);
        F = dot(Xu,Xv,2);
        G = dot(Xv,Xv,2);

        I(1,1) = E;
        I(1,2) = F;
        I(2,1) = F;
        I(2,2) = G;

        L = dot(Xuu,Nrm,2);
        M = dot(Xuv,Nrm,2);
        N = dot(Xvv,Nrm,2);

        II(1,1) = L;
        II(1,2) = M;
        II(2,1) = M;
        II(2,2) = N;

    end
    
    
end

end

function [x,w,t] = fejer1(n)
    % Fejer1 nodes: k=1/2,3/2,...,n-1/2; vector of weights: wf1
    N = (1:2:n-1)';
    l = length(N);
    m = n-l;
    K = (0:m-1)';
    v0 = [2*exp(1i*pi*K/n)./(1-4*K.^2);zeros(l+1,1)];
    v1 = v0(1:end-1)+conj(v0(end:-1:2));
    w = ifft(v1);

    k = (0:n-1)';
    k = (1/2+k);
    x = -cos(k*pi/n);
    
    t = k*pi/n;
 end