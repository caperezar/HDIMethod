function HSf = HS(f,bdy)
% Hypersingular operator computed using HOSS

    nElem = numel(bdy);

    HSf(nElem) = Density;

    for a=1:nElem

        HSf(a) = HSf(a).cctor(bdy(a));

        for b=1:nElem

            if a==b
                v = HSself(f(b),bdy(b));
            else
                v = HSfar(f(b),bdy(a),bdy(b));
            end

        end

        HSf(a).vec = HSf(a).vec+v.vec;

    end

end

function HSf = HSself(f,bdy)

    HSf = Density;
    HSf = HSf.cctor(bdy,f.dim);
    
    dfu = f.diff(bdy,'U');
    dfv = f.diff(bdy,'V');    
    dfuu = dfu.diff(bdy,'U');
    dfuv = dfv.diff(bdy,'U');
    dfvv = dfv.diff(bdy,'V');
  
    for i = 1:bdy.nPat   

        ni = bdy.nPts(i);

        for p=1:ni

            for q=1:ni

                x = bdy.Xp(q+(p-1)*ni,:,i);
                
                nrm_x = bdy.Nrm(q+(p-1)*ni,:,i);
                
                C = bdy.invA(p,q,i); 

                indX = bdy.loc_to_glob(p,q,i);

                c = zeros(size(C,1),f.dim);

                for d=1:f.dim
                
                    f0 = f.vec(indX,:);
                    fU = dfu.vec(indX,:);
                    fV = dfv.vec(indX,:);
                    fUU = dfuu.vec(indX,:);
                    fUV = dfuv.vec(indX,:);
                    fVV = dfvv.vec(indX,:);

                    rhs1 = [f0;fU;fV]; 
                    rhs2 = 0*[f0;fU;fV];
                    rhs3 = [fUU;fVV;fUV]; 
            
                    c(:,d)= C*[rhs1;rhs2;rhs3];

                end

                for j=1:bdy.nPat

                    R1 = x(1)-bdy.Xp(:,1,j);
                    R2 = x(2)-bdy.Xp(:,2,j);
                    R3 = x(3)-bdy.Xp(:,3,j);

                    HarmonicPol = [ones(size(R1)), -R1, -R2, -R3, R1.*R2, R1.*R3, R2.*R3,...
                           R1.^2 + R2.^2 - 2*R3.^2, R1.^2 - R2.^2];
                       
                    nrm_y = bdy.Nrm(:,:,j);
                    
                    R = sqrt(R1.^2+R2.^2+R3.^2)';
                    
                    invR3 = 1./R.^3/(4*pi);
                    invR5 = invR3./R.^2;
                    
                    if i==j                       
                        invR5(q+(p-1)*ni)=0.0;
                        invR3(q+(p-1)*ni)=0.0;
                    end                    
                                        
                    R_nrm_x = (nrm_x(1)*R1+nrm_x(2)*R2+nrm_x(3)*R3)';
                    R_nrm_y = dot(nrm_y,[R1,R2,R3],2)';

                    nrm_x_nrm_y = (nrm_x(1)*nrm_y(:,1)+nrm_x(2)*nrm_y(:,2)+nrm_x(3)*nrm_y(:,3))';
                                                           
                    d2Gdn2 = -3*R_nrm_x.*R_nrm_y.*invR5 + nrm_x_nrm_y.*invR3;
                    dGdn =  -R_nrm_x.*invR3;
                                        
                    Jcb = bdy.Jp(:,j).*bdy.Wp;

                    for d=1:f.dim

                        F = f.density_to_vector(j,d);
                        
                        P  = HarmonicPol*c(:,d);
                        Q  = dot(bdy.Nrm(:,:,j), [c(2,d)-c(5,d)*R2-c(6,d)*R3-2*c(8,d)*R1-2*c(9,d)*R1,...
                                          c(3,d)-c(5,d)*R1-c(7,d)*R3-2*c(8,d)*R2+2*c(9,d)*R2,...
                                          c(4,d)-c(6,d)*R1-c(7,d)*R2+4*c(8,d)*R3],2); 
                                      
                        HSf.vec(indX,d) = HSf.vec(indX,d) + dGdn*(Q.*Jcb)+d2Gdn2*((F-P).*Jcb);

                    end

                end

            end

        end
    end    
end

function HSf = HSfar(f,bdy_evl,bdy_int)

    HSf = Density;
    HSf = HSf.cctor(bdy_evl,f.dim);

    ind = [];
    outVec = [];
    
    for i = 1:bdy_evl.nPat   

        v = zeros(bdy_evl.nPts(i).^2,1);
        
        ind_aux = [];        
        
        for j=1:bdy_int.nPat

            R1 = repmat(bdy_evl.Xp(:,1,i),1,bdy_int.nPts(j)^2)-repmat(bdy_int.Xp(:,1,j)',bdy_int.nPts(i)^2,1);
            R2 = repmat(bdy_evl.Xp(:,2,i),1,bdy_int.nPts(j)^2)-repmat(bdy_int.Xp(:,2,j)',bdy_int.nPts(i)^2,1);
            R3 = repmat(bdy_evl.Xp(:,3,i),1,bdy_int.nPts(j)^2)-repmat(bdy_int.Xp(:,3,j)',bdy_int.nPts(i)^2,1);
            
            Nx1 = repmat(bdy_evl.Nrm(:,1,i),1,bdy_int.nPts(j)^2);
            Nx2 = repmat(bdy_evl.Nrm(:,2,i),1,bdy_int.nPts(j)^2);
            Nx3 = repmat(bdy_evl.Nrm(:,3,i),1,bdy_int.nPts(j)^2);
            
            Ny1 = repmat(bdy_int.Nrm(:,1,j)',bdy_int.nPts(i)^2,1);
            Ny2 = repmat(bdy_int.Nrm(:,2,j)',bdy_int.nPts(i)^2,1);
            Ny3 = repmat(bdy_int.Nrm(:,3,j)',bdy_int.nPts(i)^2,1);
            
            R = sqrt(R1.^2+R2.^2+R3.^2);
            
            G = 1./(4*pi*R);
            dG = -G./(R.^2);
            d2G = 2*G./R.^2;
                                        
            R_nrm_x = Nx1.*R1 + Nx2.*R2 + Nx3.*R3;
            R_nrm_y = Ny1.*R1 + Ny2.*R2 + Ny3.*R3;

            nrm_x_nrm_y = Nx1.*Ny1+Nx2.*Ny2+Nx3.*Ny3;
             
            d2Gdn2 = (dG-R.*d2G).*R_nrm_x.*R_nrm_y./R.^3-dG.*nrm_x_nrm_y./R;
                            
            Jcb = bdy.Jp(:,j).*bdy.Wp;
            
            F = f.density_to_vector(j,d);

            v = v + d2Gdn2*(F.*Jcb);
            
            %%
            [minR,loc]=min(R,[],2);
            ind_evl_aux =find(minR<bdy_int.tol);
            ind_aux =[ind_aux;ind_evl_aux i loc(ind_evl_aux) j]; 
            %%  

        end
        
        ind = [ind;ind_aux];
        outVec = [outVec;v];
       
    end
    
    if ~isempty(ind)
    
        dfu = f.diff(bdy_int,'U');
        dfv = f.diff(bdy_int,'V');

        dfuu = dfu.diff(bdy_int,'U');
        dfuv = dfv.diff(bdy_int,'U');
        dfvv = dfv.diff(bdy_int,'V');                    

        Nnear = size(ind,1);
                
        for i = 1:Nnear
            
            x = bdy_evl.Xp(ind(i,1),:,ind(i,2));
            nrm_x = bdy_evl.Nrm(ind(i,1),:,ind(i,2));
            
            indX = ind(i,1)+sum(bdy_evl.nPts(1:ind(i,2)-1).^2);
            
            outVec(indX)=0.0;                        
            
            indX0 = ind(i,3);
            idX0 = ind(i,4);
            
            n0 = bdy_int.nPts(idX0);
            
            if mod(indX0,n0)==0
                q = bdy_int.nPts(idX0);
            else
                q = mod(indX0,n0);
            end
            
            p = 1+(indX0-q)/ne;
            
            C = bdy_int.invA(bdy_int.chebGrid(p),bdy_int.chebGrid(q),idX0);
                               
            f0 = f.vec(indX0,:);
            fU = dfu.vec(indX0,:);
            fV = dfv.vec(indX0,:);
            fUU = dfuu.vec(indX0,:);
            fUV = dfuv.vec(indX0,:);
            fVV = dfvv.vec(indX0,:);

            rhs1 = -[f0;fU;fV]; 
            rhs2 = 0*[f0;fU;fV];
            rhs3 = -[fUU;fVV;fUV]; 

            c= C*[rhs1;rhs2;rhs3];
                                    
            for j=1:bdy.nPat

                R1 = x(1)-bdy_int.Xp(:,1,j);
                R2 = x(2)-bdy_int.Xp(:,2,j);
                R3 = x(3)-bdy_int.Xp(:,3,j);
    
                HarmonicPol = [ones(size(R1)), -R1, -R2, -R3, R1.*R2, R1.*R3, R2.*R3,...
                           R1.^2 + R2.^2 - 2*R3.^2, R1.^2 - R2.^2];
                
                nrm_y = bdy_int.Nrm(:,:,j);

                R = sqrt(R1.^2+R2.^2+R3.^2)';

                G = 1./(4*pi*R);
                dG = -G./R;
                d2G =2.*G./R.^2;

                R_nrm_x = (nrm_x(1)*R1+nrm_x(2)*R2+nrm_x(3)*R3)';
                R_nrm_y = dot(nrm_y,[R1,R2,R3],2)';

                nrm_x_nrm_y = (nrm_x(1)*nrm_y(:,1)+nrm_x(2)*nrm_y(:,2)+nrm_x(3)*nrm_y(:,3))';

                d2Gdn2 = (dG-R.*d2G).*R_nrm_x.*R_nrm_y./R.^3-dG.*nrm_x_nrm_y./R;
                dGdn = dG.*R_nrm_x./R;

                Jcb = bdy_int.Jp(:,j).*bdy_int.Wp;

                F = f.density_to_vector(j);

                P  = HarmonicPol*c;
                Q  = dot(bdy_int.Nrm(:,:,j), [c(2)-c(5)*R2-c(6)*R3-2*c(8)*R1-2*c(9)*R1,...
                                  c(3)-c(5)*R1-c(7)*R3-2*c(8)*R2+2*c(9)*R2,...
                                  c(4)-c(6)*R1-c(7)*R2+4*c(8)*R3],2); 

                outVec(indX) = outVec(indX) + dGdn*(Q.*Jcb)+d2Gdn2*((F-P).*Jcb);
                
            end
        end

    end
     
    HSf.vec = outVec;
end