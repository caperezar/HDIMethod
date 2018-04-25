function SLf = SL(f,bdy)
% Single layer operator computed using HOSS

    nElem = numel(bdy);

    SLf(nElem) = Density;

    for a=1:nElem

        SLf(a) = SLf(a).cctor(bdy(a));

        for b=1:nElem

            if a==b            
                v = SLself(f(b),bdy(b));                                
            else
                v = SLfar(f(b),bdy(a),bdy(b));
            end
            
            SLf(a).vec = SLf(a).vec+v.vec;
            
        end
        
    end

end

function SLf = SLself(f,bdy)

    SLf = Density;
    SLf = SLf.cctor(bdy,f.dim);
    
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

                    rhs1 = 0*[f0;fU;fV]; 
                    rhs2 = [f0;fU;fV];
                    rhs3 = 0*[fUU;fVV;fUV]; 

                    c(:,d)= C*[rhs1;rhs2;rhs3];

                end

                for j=1:bdy.nPat

                    R1 = x(1)-bdy.Xp(:,1,j);
                    R2 = x(2)-bdy.Xp(:,2,j);
                    R3 = x(3)-bdy.Xp(:,3,j);

                    R = sqrt(R1.^2+R2.^2+R3.^2)';
                    
                    HarmonicPol = [ones(size(R1)), -R1, -R2, -R3, R1.*R2, R1.*R3, R2.*R3,...
                               R1.^2 + R2.^2 - 2*R3.^2, R1.^2 - R2.^2];
                                      
                    G = 1./(4*pi*R);
                    dGdn = dot(bdy.Nrm(:,:,j),[R1,R2,R3],2)'.*G./(R.^2);

                    if i==j 
                        G(q+(p-1)*ni)=0.0;
                        dGdn(q+(p-1)*ni)=0.0;                       
                    end
                    
                    Jcb = bdy.Jp(:,j).*bdy.Wp;

                    for d=1:f.dim
                        
                        F = f.density_to_vector(j,d);
                        
                        P  = HarmonicPol*c(:,d);
                        Q  = dot(bdy.Nrm(:,:,j), [c(2,d)-c(5,d)*R2-c(6,d)*R3-2*c(8,d)*R1-2*c(9,d)*R1,...
                                          c(3,d)-c(5,d)*R1-c(7,d)*R3-2*c(8,d)*R2+2*c(9,d)*R2,...
                                          c(4,d)-c(6,d)*R1-c(7,d)*R2+4*c(8,d)*R3],2); 
                                      
                        SLf.vec(indX,d) = SLf.vec(indX,d) + G*((F-Q).*Jcb)+dGdn*(P.*Jcb);
                        
                    end

                end

            end

        end
    end

end

function SLf = SLfar(f,bdy_evl,bdy_int)
    SLf = Density;
    SLf = SLf.cctor(bdy_evl,1);   
    xEvl = [];
    
    for p = 1:bdy_evl.nPat
        
        xEvl = [xEvl;bdy_evl.Xp(:,:,p)];
        
    end
    
    SLf.vec = evalSL(xEvl,f,bdy_int);
           
end



