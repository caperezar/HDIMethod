classdef Sphere
    
properties
    rad=1
    center = [0 0 0];
      
    nPat=6
    
    h = 1e-4;
end

methods
    
    
    function x = X(obj,u,v,id) % parametrization
        
        N = size(u,1);
        
        One = ones(N,1);
        
        if (id==1) 

            x = [One,u,v];

        elseif(id==2) 

            x = [-u,One,v];

        elseif (id==3) 

            x = [u,v,One];

        elseif (id==4) 

            x =[-One,-u,v];

        elseif (id==5)   

            x = [u,-One,v];

        elseif (id==6) 

            x = [-u,v,-One];

        end 

        x = obj.center+obj.rad*x./repmat(sqrt(One+u.^2+v.^2),1,3);

    end 
    
    function [dxdu,dxdv,n,jac] = dX(obj,u,v,id) % derivative of the parametrization
        
        N = size(u,1);
        
        One = ones(N,1);
        
        R = sqrt(One+u.^2+v.^2);
    
        if (id==1) 

            dxdu(:,1) = -u;
            dxdv(:,1) = -v;

            dxdu(:,2) = (One+v.^2);
            dxdv(:,2) = -u.*v;

            dxdu(:,3) = -u.*v;
            dxdv(:,3) = (One+u.^2);
                       
        elseif (id==2) 

            dxdu(:,1) = -(One+v.^2);
            dxdv(:,1) = u.*v;

            dxdu(:,2) = -u;
            dxdv(:,2) = -v;

            dxdu(:,3) = -u.*v;
            dxdv(:,3) = (One+u.^2);
            
        elseif (id==3) 

            dxdu(:,1) = (One+v.^2);
            dxdv(:,1) = -u.*v;

            dxdu(:,2) = -u.*v;
            dxdv(:,2) = (One+u.^2);

            dxdu(:,3) = -u;
            dxdv(:,3) = -v;
            
        elseif (id==4) 

            dxdu(:,1) = u;
            dxdv(:,1) = v;

            dxdu(:,2) = -(One+v.^2);
            dxdv(:,2) = u.*v;

            dxdu(:,3) = -u.*v;
            dxdv(:,3) = (One+u.^2);
            
        elseif (id==5) 

            dxdu(:,1) = (One+v.^2);
            dxdv(:,1) = -u.*v;

            dxdu(:,2) = u;
            dxdv(:,2) = v;

            dxdu(:,3) = -u.*v;
            dxdv(:,3) = (One+u.^2);
            
        elseif (id==6) 

            dxdu(:,1) = -(One+v.^2);
            dxdv(:,1) = u.*v;

            dxdu(:,2) = -u.*v;
            dxdv(:,2) = (One+u.^2);

            dxdu(:,3) = u;
            dxdv(:,3) = v;
            
            

        end 
        
        dxdu = obj.rad*dxdu./repmat(R.^3,1,3);
        dxdv = obj.rad*dxdv./repmat(R.^3,1,3);
        
        e1_e2 = cross(dxdu(:,:),dxdv(:,:));   % Normal vector    
        jac = vecnorm(e1_e2,2,2) ; % Unit normal
        n(:,:)= e1_e2 ./ jac ; % Unit normal

    end
       
    function [d2xdu2,d2xdudv,d2xdv2]=d2X(obj,u,v,id) % second derivative of the parametrization
        
        N = size(u,1);
        
        One = ones(N,1);
        
        R = sqrt(One+u.^2+v.^2);

        if (id==1) 

            d2xdu2(:,1) = (2.0*u.^2-v.^2-One);
            d2xdudv(:,1) = 3.0*u.*v;
            d2xdv2(:,1) = (2.0*v.^2-u.^2-One);

            d2xdu2(:,2) = -3.0*u.*(v.^2+One);
            d2xdudv(:,2) = v.*(2.0*u.^2-v.^2-One);
            d2xdv2(:,2) = -u.*(u.^2-2.0*v.^2+One);

            d2xdu2(:,3) = v.*(2.0*u.^2-v.^2-One);
            d2xdudv(:,3) =-u.*(u.^2-2.0*v.^2+1);
            d2xdv2(:,3) = -3.0*v.*(u.^2+One);
    
        elseif (id==2) 

            d2xdu2(:,2) = (2.0*u.^2-v.^2-One);
            d2xdudv(:,2) = 3.0*u.*v;
            d2xdv2(:,2) = (2.0*v.^2-u.^2-One);

            d2xdu2(:,1) = 3.0*u.*(v.^2+One);
            d2xdudv(:,1) = -v.*(2.0*u.^2-v.^2-One);
            d2xdv2(:,1) = u.*(u.^2-2.0*v.^2+One);

            d2xdu2(:,3) = v.*(2.0*u.^2-v.^2-One);
            d2xdudv(:,3) =-u.*(u.^2-2.0*v.^2+One);
            d2xdv2(:,3) = -3.0*v.*(u.^2+One);

        elseif (id==3) 

            d2xdu2(:,3) = (2.0*u.^2-v.^2-One);
            d2xdudv(:,3) = 3.0*u.*v;
            d2xdv2(:,3) = (2.0*v.^2-u.^2-One);

            d2xdu2(:,1) = -3.0*u.*(v.^2+One);
            d2xdudv(:,1) = v.*(2.0*u.^2-v.^2-One);
            d2xdv2(:,1) = -u.*(u.^2-2.0*v.^2+One);

            d2xdu2(:,2) = v.*(2.0*u.^2-v.^2-One);
            d2xdudv(:,2) =-u.*(u.^2-2.0*v.^2+One);
            d2xdv2(:,2) = -3.0*v.*(u.^2+One);

        elseif (id==4) 

            d2xdu2(:,1) = -(2.0*u.^2-v.^2-One);
            d2xdudv(:,1) = -3.0*u.*v;
            d2xdv2(:,1) = -(2.0*v.^2-u.^2-One);

            d2xdu2(:,2) = 3.0*u.*(v.^2+One);
            d2xdudv(:,2) = -v.*(2.0*u.^2-v.^2-One);
            d2xdv2(:,2) = u.*(u.^2-2.0*v.^2+One);

            d2xdu2(:,3) = v.*(2.0*u.^2-v.^2-One);
            d2xdudv(:,3) =-u.*(u.^2-2.0*v.^2+One);
            d2xdv2(:,3) = -3.0*v.*(u.^2+One);

        elseif (id==5)   

            d2xdu2(:,2) = -(2.0*u.^2-v.^2-One);
            d2xdudv(:,2) = -3.0*u.*v;
            d2xdv2(:,2) = -(2.0*v.^2-u.^2-One);

            d2xdu2(:,1) = -3.0*u.*(v.^2+One);
            d2xdudv(:,1) = v.*(2.0*u.^2-v.^2-One);
            d2xdv2(:,1) = -u.*(u.^2-2.0*v.^2+One);

            d2xdu2(:,3) = v.*(2.0*u.^2-v.^2-One);
            d2xdudv(:,3) =-u.*(u.^2-2.0*v.^2+One);
            d2xdv2(:,3) = -3.0*v.*(u.^2+One);

        elseif (id==6) 

            d2xdu2(:,3) = -(2.0*u.^2-v.^2-One);
            d2xdudv(:,3) = -3.0*u.*v;
            d2xdv2(:,3) = -(2.0*v.^2-u.^2-One);

            d2xdu2(:,1) = 3.0*u.*(v.^2+One);
            d2xdudv(:,1) = -v.*(2.0*u.^2-v.^2-One);
            d2xdv2(:,1) = u.*(u.^2-2.0*v.^2+One);

            d2xdu2(:,2) = v.*(2.0*u.^2-v.^2-One);
            d2xdudv(:,2) =-u.*(u.^2-2.0*v.^2+One);
            d2xdv2(:,2) = -3.0*v.*(u.^2+One);

        end

        d2xdu2 = obj.rad*d2xdu2./repmat(R.^5,1,3);
        d2xdv2 = obj.rad*d2xdv2./repmat(R.^5,1,3);
        d2xdudv = obj.rad*d2xdudv./repmat(R.^5,1,3);
    
    end 
    
    
    function [nU,nV] = dn(obj,u,v,id) % derivative of the parametrization
%        h = 1e-7;        
%        [~,~,n] = dX(obj,u,v,id);
%        [~,~,npu] = dX(obj,u+h,v,id);
%        [~,~,npv] = dX(obj,u,v+h,id);
%        n0U = (npu-n)/h;
%        n0V = (npv-n)/h;       
       
             
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
        
%         norm(n0U-nU)/norm(nU)
%         norm(n0V-nV)/norm(nV)
       
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

