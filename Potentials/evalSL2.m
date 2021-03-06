function pot = evalSL2(xEvl,f,bdy)

Np = size(xEvl,1);
pot = zeros(Np,1);

Nobs = numel(bdy);

for o=1:Nobs
    
    for j = 1:bdy(o).nPat

        nj = bdy(o).nPts(j);

        R1 = repmat(xEvl(:,1),1,nj^2)-repmat(bdy(o).Xp(:,1,j)',Np,1);
        R2 = repmat(xEvl(:,2),1,nj^2)-repmat(bdy(o).Xp(:,2,j)',Np,1);
        R3 = repmat(xEvl(:,3),1,nj^2)-repmat(bdy(o).Xp(:,3,j)',Np,1);   

        R = sqrt(R1.^2+R2.^2+R3.^2);

        G = 1./(4*pi*R);

        Jcb = bdy(o).Jp(:,j).*bdy(o).Wp;

        F = f(o).density_to_vector(j);

        pot = pot + G*(F.*Jcb);
             
    end    
end


end