classdef Cushion
    
properties
    rad=1
    center = [0 0 0];
      
    nPat=6
    
    h = 1e-3;
    rot = [0 pi/4 pi/3];
%     rot = [0 0 0];
    
end

methods
    
    
    function x = X(obj,u,v,id) % parametrization
        Rx = [1 0 0;0 cos(obj.rot(1)) sin(obj.rot(1));0 -sin(obj.rot(1)) cos(obj.rot(1))];
        Ry = [cos(obj.rot(2)) 0 -sin(obj.rot(2));0 1 0;sin(obj.rot(2)) 0 cos(obj.rot(2))];
        Rz= [cos(obj.rot(3)) sin(obj.rot(3)) 0;-sin(obj.rot(3)) cos(obj.rot(3)) 0;0 0 1];
        R = Rz*Ry*Rx;
        
        sph = Sphere;
    
        y = sph.X(u,v,id);
        
        y = y/sph.rad-sph.center;
        
        [th,phi] = cart2sph(y(:,1),y(:,2),y(:,3));

        r = sqrt(0.8+0.5*(cos(2*th)-1).*(cos(4*phi)-1));

        x(:,1) = r.*cos(th).*cos(phi);
        x(:,2) = r.*sin(th).*cos(phi);
        x(:,3) = r.*sin(phi);
        x= x*R';
        x  = obj.rad*x+obj.center;
                  
    end 
    
    function [dxdu,dxdv,n,jac] = dX(obj,u,v,id) % derivative of the parametrization
%         h0=1e-7;
        h0 = obj.h;

        UL = obj.X(u-h0,v,id);
        UR = obj.X(u+h0,v,id);
        
        dxdu = (UR-UL)/(2*h0);
        
        UL = obj.X(u,v-h0,id);
        UR = obj.X(u,v+h0,id);
        
        dxdv = (UR-UL)/(2*h0);
        
%         dxdu = obj.rad*dxdu;
%         dxdv = obj.rad*dxdv;
        
        e1_e2 = cross(dxdu(:,:),dxdv(:,:));   % Normal vector  
        jac  = vecnorm(e1_e2,2,2);
        n(:,:)= e1_e2 ./ jac ; % Unit normal

    end
       
    function [d2xdu2,d2xdudv,d2xdv2]=d2X(obj,u,v,id) % second derivative of the parametrization
        h0 = obj.h;
        

        [UL,VL] = obj.dX(u-h0,v,id);
        [UR,VR] = obj.dX(u+h0,v,id);

        d2xdu2 = (UR-UL)/(2*h0);
        d2xdudv = (VR-VL)/(2*h0);

        [UL,VL]= obj.dX(u,v-h0,id);
        [UR,VR]= obj.dX(u,v+h0,id);

        d2xdv2 = (VR-VL)/(2*h0);
        d2xdudv = 0.5*d2xdudv+0.5*(UR-UL)/(2*h0);
    
    
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
    
    function [u,v,isin]=invX(obj,x,id) % inverse of the parametrization
        
        x = x-obj.center;
        
        if (id==1) 

            u = x(:,2)./x(:,1);
            v = x(:,3)./x(:,1);

        elseif (id==2) 

            u = x(:,1)./x(:,2);
            v = x(:,3)./x(:,2);

        elseif (id==3) 

            u = x(:,1)./x(:,3);
            v = x(:,2)./x(:,3);

        elseif (id==4) 

            u = -x(:,2)./x(:,1);
            v = -x(:,3)./x(:,1);

        elseif (id==5)

            u = -x(:,1)./x(:,2);
            v = -x(:,3)./x(:,2);

        elseif (id==6) 

            u = -x(:,1)./x(:,3);
            v = -x(:,2)./x(:,3);

        end

        isin = and(abs(u)<=1,abs(v)<=1); 

    end 
    


end

end