classdef Plane
    
properties
    dir = [1 0 0;0 1 0];
    center = [0 0 0];
    L = [2 2];
    N = [20 20];   
    Xp
end

methods
        
    function obj = cctor(obj,cntr,Dir,L,N) 
        % cntr:  center of the plane
        % Dir : director vectors of the plane (rows on the matrix Dir)
        % L   : length in each direction
        % N   : number of discretization points in each direction 
        
        obj.dir = Dir;
        obj.center=cntr;
        obj.L=L;
        obj.N=N;
        
        u = linspace(-L(1)/2,L(1)/2,N(1));
        v = linspace(-L(2)/2,L(2)/2,N(2));
        
        [U,V] = meshgrid(u,v);
        
        uVec=reshape(U,N(1)*N(2),1);
        vVec=reshape(V,N(2)*N(1),1);
        
        obj.Xp(:,1) = cntr(1) + uVec*Dir(1,1) + vVec*Dir(2,1);
        obj.Xp(:,2) = cntr(2) + uVec*Dir(1,2) + vVec*Dir(2,2);
        obj.Xp(:,3) = cntr(3) + uVec*Dir(1,3) + vVec*Dir(2,3);
                
    end
    
    function []=plot(obj,vec,nfig,opt)
        % vec: vector containing the values of the field on the plane
        % nfig: figure number where the plane will be displayed
        % opt: function of the field that will to be displayed
        n1 = obj.N(1);
        n2 = obj.N(2);
        X=reshape(obj.Xp(:,1),n2,n1);
        Y=reshape(obj.Xp(:,2),n2,n1);
        Z=reshape(obj.Xp(:,3),n2,n1);
        V = reshape(vec,n2,n1);
        figure(nfig);hold on
        if strcmp(opt,'real')
            surf(X,Y,Z,real(V),'edgecolor','none');
            colormap(brewermap([],'*RdYlBu'))
        elseif strcmp(opt,'imag')    
            surf(X,Y,Z,imag(V),'edgecolor','none');
            colormap(brewermap([],'*RdYlBu'))
        elseif strcmp(opt,'abs')    
            surf(X,Y,Z,abs(V),'edgecolor','none');
            colormap(brewermap([],'*YlGnBu'))
        elseif strcmp(opt,'log10')    
            surf(X,Y,Z,log10(abs(V)),'edgecolor','none');
%             colormap(brewermap([],'*YlGnBu')) 
            colormap hot
        end            
    end
    


end

end
