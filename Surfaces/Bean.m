classdef Bean
    
properties
    
    rad=1
    center = [0 0 0];
    
    nPat=6
        
    h = 1e-3;
end

methods
    
    function x = X(obj,u,v,id) % parametrization
        
        sph = Sphere;
        
        alpha1 = 0.3;
        alpha2 = 0.4;
        alpha3 = 0.1;

        a = 0.8;
        b = 0.8;
    
        y = sph.X(u,v,id);
        
        y = y/sph.rad-sph.center;
    
        x(:,1) = a*sqrt(1.0-alpha3*cos(pi*y(:,3))).*y(:,1);
        x(:,2) =-alpha1*cos(pi*y(:,3))+b*sqrt(1.0-alpha2*cos(pi*y(:,3))).*y(:,2);
        x(:,3) = y(:,3);
        
        x  = obj.rad*x+obj.center;
    
    end 
    
    function [dxdu,dxdv,n,jac] = dX(obj,u,v,id) % derivative of the parametrization
        
        sph = Sphere;  
        
        y =sph.X(u,v,id);
        
        [dydu,dydv] = sph.dX(u,v,id);

        y = y/sph.rad-sph.center;

        dydu = dydu/sph.rad;

        dydv = dydv/sph.rad;

        alpha1 = 0.3;
        alpha2 = 0.4;
        alpha3 = 0.1;
        a = 0.8;
        b = 0.8;
        ONE = ones(size(u,1),1);
        ZERO = zeros(size(u,1),1);
            
        dx1dy(:,1) = a*sqrt(1.0-alpha3*cos(pi*y(:,3)));
        dx1dy(:,2) = ZERO;
        dx1dy(:,3) = 0.5*a*alpha3*pi*y(:,1).*sin(pi*y(:,3))./sqrt(1.0-alpha3*cos(pi*y(:,3)));

        dx2dy(:,1) = ZERO;
        dx2dy(:,2) = b*sqrt(1.0-alpha2*cos(pi*y(:,3)));
        dx2dy(:,3) = alpha1*pi*sin(pi*y(:,3))+0.5*b*alpha2*pi*y(:,2).*sin(pi*y(:,3))./sqrt(1.0-alpha2*cos(pi*y(:,3)));
     

        dx3dy(:,1) = ZERO;
        dx3dy(:,2) = ZERO;
        dx3dy(:,3) = ONE;
        
        dxdu(:,1) = dot(dx1dy,dydu,2);
        dxdu(:,2) = dot(dx2dy,dydu,2);
        dxdu(:,3) = dot(dx3dy,dydu,2);

        dxdv(:,1) = dot(dx1dy,dydv,2);
        dxdv(:,2) = dot(dx2dy,dydv,2);
        dxdv(:,3) = dot(dx3dy,dydv,2);
        
        
        dxdu = obj.rad*dxdu;
        dxdv = obj.rad*dxdv;
        
        e1_e2 = cross(dxdu(:,:),dxdv(:,:));   % Normal vector  
        jac  = vecnorm(e1_e2,2,2);
        n(:,:)= e1_e2 ./ jac ; % Unit normal

    end
       
    function [d2xdu2,d2xdudv,d2xdv2]=d2X(obj,u,v,id) % second derivative of the parametrization
        
        h = 1.0e-5;

        [UL,VL] = obj.dX(u-h,v,id);
        [UR,VR] = obj.dX(u+h,v,id);

        d2xdu2 = (UR-UL)/(2*h);
        d2xdudv = (VR-VL)/(2*h);

        [UL,VL]= obj.dX(u,v-h,id);
        [UR,VR]= obj.dX(u,v+h,id);

        d2xdv2 = (VR-VL)/(2*h);
        d2xdudv = 0.5*d2xdudv+0.5*(UR-UL)/(2*h);
    
    end 
    
    
    function [nU,nV] = dn(obj,u,v,id) % derivative of the parametrization
%         h = 1e-7;        
%        [~,~,n] = dX(obj,u,v,id);
%        [~,~,npu] = dX(obj,u+h,v,id);
%        [~,~,npv] = dX(obj,u,v+h,id);
%        nU = (npu-n)/h;
%        nV = (npv-n)/h;       

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

