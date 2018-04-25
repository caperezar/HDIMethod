classdef Cube
    
properties
    
    rad=1
    center = [0 0 0];
    
    nPat=6
    h = 1e-4;
    
    stretch = [1 0.25 1];
    rot = [0 -pi/3 pi/7];
end

methods
    
    
    function x = X(obj,u,v,id) % parametrization
        
        
        
        ONE = ones(size(u,1),1);
        
        if (id==1) 

            x = [ONE,u,v];

        elseif (id==2) 

            x = [-u,ONE,v];

        elseif (id==3) 

            x = [u,v,ONE];
    
        elseif (id==4) 

            x = [-ONE,-u,v];

        elseif (id==5)   

            x = [u,-ONE,v];

        elseif (id==6) 

            x = [-u,v,-ONE];

        end 
        
        Rx = [1 0 0;0 cos(obj.rot(1)) sin(obj.rot(1));0 -sin(obj.rot(1)) cos(obj.rot(1))];
        Ry = [cos(obj.rot(2)) 0 -sin(obj.rot(2));0 1 0;sin(obj.rot(2)) 0 cos(obj.rot(2))];
        Rz= [cos(obj.rot(3)) sin(obj.rot(3)) 0;-sin(obj.rot(3)) cos(obj.rot(3)) 0;0 0 1];
        Rrot = Rz*Ry*Rx;
        
        R = diag(obj.stretch);
                
        x = x*(Rrot'*R*Rrot)';

    end 
    
    function [dxdu,dxdv,n,jac] = dX(obj,u,v,id) % derivative of the parametrization
        
        ZERO = zeros(size(u,1),1);
        ONE = ones(size(u,1),1);
        
        if (id==1) 

            dxdu = [ZERO,ONE,ZERO];
            dxdv = [ZERO,ZERO,ONE];
            
        elseif (id==2) 

            dxdu = [-ONE,ZERO,ZERO];
            dxdv = [ZERO,ZERO,ONE];
            
        elseif (id==3) 

            dxdu = [ONE,ZERO,ZERO];
            dxdv = [ZERO,ONE,ZERO];
            
        elseif (id==4) 

            dxdu = [ZERO,-ONE,ZERO];
            dxdv = [ZERO,ZERO,ONE];
            
        elseif (id==5)   

            dxdu = [ONE,ZERO,ZERO];
            dxdv = [ZERO,ZERO,ONE];
            
        elseif (id==6) 

            dxdu = [-ONE,ZERO,ZERO];
            dxdv = [ZERO,ONE,ZERO];

        end 
        
        Rx = [1 0 0;0 cos(obj.rot(1)) sin(obj.rot(1));0 -sin(obj.rot(1)) cos(obj.rot(1))];
        Ry = [cos(obj.rot(2)) 0 -sin(obj.rot(2));0 1 0;sin(obj.rot(2)) 0 cos(obj.rot(2))];
        Rz= [cos(obj.rot(3)) sin(obj.rot(3)) 0;-sin(obj.rot(3)) cos(obj.rot(3)) 0;0 0 1];
        Rrot = Rz*Ry*Rx;
        
        R = diag(obj.stretch);
        
        
        dxdu = dxdu*(Rrot'*R*Rrot)';
        dxdv = dxdv*(Rrot'*R*Rrot)';
        
        e1_e2 = cross(dxdu(:,:),dxdv(:,:));   % Normal vector    
        jac = vecnorm(e1_e2,2,2);
        n(:,:)= e1_e2 ./jac; % Unit normal

    end
       
    function [d2xdu2,d2xdudv,d2xdv2]=d2X(obj,u,v,id) % second derivative of the parametrization
                
        d2xdu2 = zeros(size(u,1),3);
        d2xdudv = zeros(size(u,1),3);
        d2xdv2 = zeros(size(u,1),3);
                
    end 
    
    
    function [nU,nV] = dn(obj,u,v,id) % derivative of the parametrization
        nU = zeros(size(u,1),3);
        nV= zeros(size(u,1),3);
        
    end
    
end

end
