function Ac = my_cobec( Y,C )
%COMSSPACE Summary of this function goes here
%   Detailed explanation goes here
% defopts=struct('c',1,'maxiter',500,'ctol',1e-3,'PCAdim',[]);
maxiter = 500;
ctol = 1e-3;

NRows=size(Y,1);
N=size(Y,3);

Ac=ccak_init(Y(:,:,1),Y(:,:,2),C);

%% iterations
x=cell(1,N);
for it=1:maxiter
    c0=Ac;
    
    c2=zeros(NRows,C);
    for n=1:N
        x{n}=Y(:,:,n)'*Ac;
        c2=c2+Y(:,:,n)*x{n};
    end
   
    %% SVDS
    Ac=c2*(c2'*c2)^-.5;
    
    %% stop
    if mean(abs(diag(Ac'*c0)))>1-ctol
        break;
    end
end
end

