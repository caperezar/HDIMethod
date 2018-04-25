function [densOut,it] = GMRESsolve(densIn,bdy)
% GMRES parameters
tol = 10^(-8);
maxit = 50;

Nobs = numel(bdy);

fvec = [];
for o=1:Nobs
    fvec = [fvec;densIn(o).vec];    
end 


[densVec,~,~,it,~] = gmres(@(vec)Map(vec,densIn,bdy),fvec,[],tol,maxit);

densOut(1:Nobs) = Density; 

for o=1:Nobs
    densOut(o) = densOut(o).cctor(bdy(o));
    densOut(o).vec = densVec(1:bdy(o).nPtot);
    densVec(1:bdy(o).nPtot)=[];
end

end

% Apply operator:
function vecOut = Map(vecIn,dens,bdy)
    Nobs  = numel(bdy);
    for o = 1:Nobs
        dens(o).vec = vecIn(1:bdy(o).nPtot);
        vecIn(1:bdy(o).nPtot) = [];
    end
    
    dlf = DL(dens,bdy);
    
    vecOut = [];
    for o = 1:Nobs
        vecOut = [vecOut;dlf(o).vec-dens(o).vec/2];        
    end
    
end
