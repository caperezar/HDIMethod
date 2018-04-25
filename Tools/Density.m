classdef Density
    
properties    
     nPat
     vec
     vecSize
     nPts
     dim
end

methods
    
    function obj=cctor(obj,bdy,dim)
        if nargin==2
            dim=1;
        end
        obj.dim = dim;
            
        obj.nPts = bdy.nPts;
        
        obj.nPat = bdy.nPat;
        
        obj.vecSize = sum(bdy.nPts.^2);
        
        obj.vec = zeros(obj.vecSize,dim);
        
    end
    
    function v=density_to_vector(obj,id,dim)
        if nargin==2
            dim=1;
        end
        
        N0=0;
        
        if (id>1) 
        
            N0 = sum(obj.nPts(1:id-1).^2);

        end

        n = obj.nPts(id);
        
        v = obj.vec(N0+1:N0+n^2,dim);
        
    end 
    
    
    function mat=density_to_matrix(obj,id,dim)
        if nargin==2
            dim=1;
        end
        
        N0=0;
        
        if (id>1) 
        
            N0 = sum(obj.nPts(1:id-1).^2);

        end

        n = obj.nPts(id);
        
        
        v = obj.vec(N0+1:N0+n^2,dim);
        
        mat = reshape(v,n,n);

    end 
   
    
    function obj = matrix_to_density(obj,id,mat,dim)
        if nargin==3
            dim=1;
        end
        n = obj.nPts(id);

        N0 = 0;
        
        if (id>1) 
             N0 = sum(obj.nPts(1:id-1).^2);
        end

        v = reshape(mat,n^2,1);

        obj.vec(N0+1:N0+n^2,dim)=v;   
        
    end
    
    function obj = vector_to_density(obj,id,v,dim)
        
        if nargin==3
            dim=1;
        end
        
        n = obj.nPts(id);

        N0 = 0;
        
        if (id>1) 
             N0 = sum(obj.nPts(1:id-1).^2);
        end

        obj.vec(N0+1:N0+n^2,dim)=v;   
        
    end
    
    function []=plot(obj,bdy,nfig,opt,dim)
         if nargin==4
            dim=1;
            
         end
         if nargin==3
             opt='real';
         end
        
        figure(nfig);hold on;        
        for j=1:bdy.nPat
            n = bdy.nPts(j);
            X=reshape(bdy.Xp(:,1,j),n,n);
            Y=reshape(bdy.Xp(:,2,j),n,n);
            Z=reshape(bdy.Xp(:,3,j),n,n);
            V = obj.density_to_matrix(j,dim);
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
                colormap hot
%                 colormap(brewermap([],'*YlGnBu'));                
            end
            
        end
        axis tight
        axis equal
        hold off        
%         light; lighting phong;% camzoom(1.3);
%         caxis([min(min(realobj.vec)) max(max(obj.vec))])
%         camlight('left');
%         camlight('right');
%         camlight('headlight','infinite');
%         camlight('left','infinite');
%         camlight('right','infinite');
        
        shading interp
    end
    
    function gradf = grad(obj,bdy,basis)
        % This function computes the surface gradient of the density f
        % (Nedelec's book formula (2.5.176))

        df1 = obj.diff(bdy,'U');
        df2 = obj.diff(bdy,'V');

        % basis to the tangent plane
        e1=[];e2=[];
        for p=1:bdy.nPat
            e1 = [e1;bdy.dXdu(:,:,p)];
            e2 = [e2;bdy.dXdv(:,:,p)];            
        end

        % Coefficients of the metric tensor
        g11 = dot(e1,e1,2); 
        g12 = dot(e1,e2,2);
        g22 = dot(e2,e2,2);

        det_g = g11.*g22-g12.^2;

        % Contravariant basis
        ep1 = (e1.*g22-e2.*g12)./det_g; 
        ep2 = (-e1.*g12+e2.*g11)./det_g;

        % Surface gradient
        if nargin ==2 || strcmp(basis,'cartesian')
            gradf = Density;
            gradf = gradf.cctor(bdy,3);
            gradf.vec = ep1.*repmat(df1.vec,1,3)+ep2.*repmat(df2.vec,1,3);

        elseif strcmp(basis,'covariant')    
            gradf = Density;
            gradf = gradf.cctor(bdy,2);
            gradf.vec =[df1.vec./det_g,df2.vec/deg_g];
        end
    end
    
    function rf = curlScalar(obj,bdy)
    % This function computes the curl scalar of the vector density v which is
    % defined as curl v dot normal (Nedelec's book formula (2.5.206))

    % Note: v is assumed to be given in the Cartesian basis

        % natural basis
        e1=[];e2=[];
        for p=1:bdy.nPat
            e1 = [e1;bdy.dXdu(:,:,p)];
            e2 = [e2;bdy.dXdv(:,:,p)];            
        end

        % Coefficients of the metric tensor
        g11 = dot(e1,e1,2); 
        g12 = dot(e1,e2,2);
        g22 = dot(e2,e2,2);

        det_g = g11.*g22-g12.^2;

        % Components of the vector field in the covariant basis
        v1 =Density;v1=v1.cctor(bdy,1);
        v2 =Density;v2=v2.cctor(bdy,1);

        % components in the covariant basis
        v1.vec = dot(e1,obj.vec,2);
        v2.vec = dot(e2,obj.vec,2);

        dv1 = v2.diff(bdy,'U');
        dv2 = v1.diff(bdy,'V');   

        %  
        rf = Density;
        rf = rf.cctor(bdy,1);

        rf.vec = (dv1.vec-dv2.vec)./sqrt(det_g);
        % rf.vec = 1./sqrt(det_g);


    end

    
    function CVf = curlVector(obj,bdy,basis)
    % This function computes the curl vector of the density f which is
    % define as grad_surf f x normal (Nedelec's book formula (2.5.177))
        df1 = obj.diff(bdy,'V');
        df2 = obj.diff(bdy,'U');

        % basis to the tangent plane
        e1=[];e2=[];
        for p=1:bdy.nPat
            e1 = [e1;bdy.dXdu(:,:,p)];
            e2 = [e2;bdy.dXdv(:,:,p)];            
        end

        % Coefficients of the metric tensor
        g11 = dot(e1,e1,2); 
        g12 = dot(e1,e2,2);
        g22 = dot(e2,e2,2);

        sq_det_g = sqrt(g11.*g22-g12.^2);

        if nargin ==2 || strcmp(basis,'cartesian')

            CVf = Density;
            CVf = CVf.cctor(bdy,3);

            CVf.vec = e1.*repmat(df1.vec./sq_det_g,1,3)-e2.*repmat(df2.vec./sq_det_g,1,3);

        elseif strcmp(basis,'natural')    

            CVf = Density;
            CVf = CVf.cctor(bdy,2);
            CVf.vec = [df1.vec./sq_det_g,-df2.vec./sq_det_g];

        end
    end
    
    function divf = div(obj,bdy)
    % This function computes the surface divergence of the vector density v
    % (Nedelec's book formula (2.5.205))

    % basis to the tangent plane
        e1=[];e2=[];
        for p=1:bdy.nPat
            e1 = [e1;bdy.dXdu(:,:,p)];
            e2 = [e2;bdy.dXdv(:,:,p)];            
        end

        % Coefficients of the metric tensor
        g11 = dot(e1,e1,2); 
        g12 = dot(e1,e2,2);
        g22 = dot(e2,e2,2);

        det_g = g11.*g22-g12.^2;
        sq_det_g = sqrt(det_g);

        % Contravariant basis
        ep1 = (e1.*g22-e2.*g12)./det_g; 
        ep2 = (-e1.*g12+e2.*g11)./det_g;


        v1 =Density;v1 =v1.cctor(bdy,1);
        v1.vec = dot(obj.vec,ep1,2).*sq_det_g;
        v1 = v1.diff(bdy,'U');

        v2 =Density;v2 =v2.cctor(bdy,1);
        v2.vec = dot(obj.vec,ep2,2).*sq_det_g;
        v2 = v2.diff(bdy,'V');


        divf = Density;
        divf = divf.cctor(bdy,1);

        divf.vec = (v1.vec + v2.vec)./sq_det_g; 

    end
    
    function dg = diff(obj,bdy,dir)

      d = obj.dim;

      L  = pi;

      dg = Density;

      dg = dg.cctor(bdy,d);

      for d=1:d

          for p  = 1:obj.nPat

             np = bdy.nPts(p);

             A = obj.density_to_matrix(p,d);

             if strcmp(dir,'U') 

                dAdU = FFTDiff1DEven(A.',L);

                dAdU = dAdU./repmat(sin(bdy.tGrid),1,np);

                dg = dg.matrix_to_density(p,dAdU.',d);

             elseif strcmp(dir,'V') 

                dAdV = FFTDiff1DEven(A,L);

                dAdV = dAdV./repmat(sin(bdy.tGrid),1,np);

                dg = dg.matrix_to_density(p,dAdV,d);

             end

          end
      end

    end
    
    function I=integral(obj,bdy,dim)
        jp=[];
        for p=1:bdy.nPat    
            jp = [jp;bdy.Jp(:,p).*bdy.Wp];
        end
        I= sum(jp.*obj.vec(:,dim));
    end
       
end

end

function dx = FFTDiff1DEven(x,L)
  
  N = size(x,1);
  Nvecs = size(x,2);
  Nx   = 2*N;

  x0_in=zeros(Nx,Nvecs);
  
  k = zeros(Nx,Nvecs);

  
  x0_in((N+1):2*N,:) = x;
  x0_in(1:N,:) = x(N:-1:1,:);

  x0 = fft(x0_in);
  
  k(1:Nx/2,:)    = repmat((0:Nx/2-1)',1,Nvecs);
  k(Nx/2+1,:)    = zeros(1,Nvecs);
  k(Nx/2+2:Nx,:) = repmat((-Nx/2+1:-1)',1,Nvecs);

  dx0 = 2i*pi/L*x0.*k;

  dy0 = ifft(dx0);
  
  dx = 2*N*dy0((N+1):2*N,:)/(2*Nx);

end

