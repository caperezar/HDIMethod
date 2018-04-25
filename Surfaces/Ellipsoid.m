classdef Ellipsoid
    
properties
    
    rad=1
    center = [0 0 0];
    a=1.0;
    b=0.375;
    c=0.5;   
    nPat=6;
    h = 1e-4;
end

methods
    
  
    
    function x = X(obj,u,v,id) % parametrization
        
        sph = Sphere;
        
        y = sph.X(u,v,id);
        
        y = y/sph.rad-sph.center;
    
        x = [obj.a*y(:,1),obj.b*y(:,2),obj.c*y(:,3)];
        
        x = obj.rad*x+obj.center;
    
    end 
    
    function [dxdu,dxdv,n,jac] = dX(obj,u,v,id) % derivative of the parametrization
        
        sph = Sphere;  
        
        [dydu,dydv] = sph.dX(u,v,id);
        
        dydu = dydu/sph.rad;

        dydv = dydv/sph.rad;
        
        dxdu = [obj.a*dydu(:,1),obj.b*dydu(:,2),obj.c*dydu(:,3)];
        dxdv = [obj.a*dydv(:,1),obj.b*dydv(:,2),obj.c*dydv(:,3)];
        
        dxdu=obj.rad*dxdu;
        dxdv=obj.rad*dxdv;
        
        e1_e2 = cross(dxdu(:,:),dxdv(:,:));   % Normal vector    
        jac = vecnorm(e1_e2,2,2);
        n(:,:)= e1_e2 ./ jac; % Unit normal

    end
       
    function [d2xdu2,d2xdudv,d2xdv2]=d2X(obj,u,v,id) % second derivative of the parametrization
        sph= Sphere;
        [d2ydu2,d2ydudv,d2ydv2] = sph.d2X(u,v,id);
        
        d2xdu2 = obj.rad*[obj.a*d2ydu2(:,1),obj.b*d2ydu2(:,2),obj.c*d2ydu2(:,3)];
        d2xdv2 = obj.rad*[obj.a*d2ydv2(:,1),obj.b*d2ydv2(:,2),obj.c*d2ydv2(:,3)];
        d2xdudv = obj.rad*[obj.a*d2ydudv(:,1),obj.b*d2ydudv(:,2),obj.c*d2ydudv(:,3)];
        
    end 
    
    
    function [nU,nV] = dn(obj,u,v,id) % derivative of the parametrization
%         h = 1e-7;        
%        [~,~,n] = dX(obj,u,v,id);
%        [~,~,npu] = dX(obj,u+h,v,id);
%        [~,~,npv] = dX(obj,u,v+h,id);
%        nU = (npu-n)/h;
%        nV = (npv-n)/h;
%        
       [xu,xv,n,jac] = obj.dX(u,v,id);
       [xuu,xuv,xvv] = obj.d2X(u,v,id);
       g11= dot(xu,xu,2);
       g12= dot(xu,xv,2);
       g22= dot(xv,xv,2);
       
       det_g = jac.^2;

        % Contravariant basis
        ep1 = (xu.*g22-xv.*g12)./det_g; 
        ep2 = (-xu.*g12+xv.*g11)./det_g;
        
        alpha_U = -dot(n,xuu,2);
        beta_U = -dot(n,xuv,2);
       
        nU = alpha_U.*ep1+beta_U.*ep2;
        
        alpha_V = beta_U;
        beta_V = -dot(n,xvv,2);
        
        nV = alpha_V.*ep1+beta_V.*ep2;
        
    end
    

   
    
end

end

