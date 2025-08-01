function [ D, Xn, DnXn, DXn ] = my_cobe( Yn ,C, Rn, n_iter) 

%% COBE Dictionary Learning Algorithm

% formulation -->        Yn = DXn + DnXn'
% input -->    Yn = a cell containing all the matrices --> size(Yn) = 1 x n (n = number of matrices) 
%              Yn{i} = a matrix of size d x m , d = dimension of features, m = number of vectors in matrix Yn{i} 
%              C = no. of atoms in common dictionary D
%              nor = 1      zscore normalization on Yn 
%                  = 0      no normalization (Default)
% output -->   D = Common dictionary --> size = m x C (m is the dimention)
[d,m] = size(Yn{1});
if nargin == 2
    Rn = m;
    n_iter = 100;
end

n_mat = length(Yn);
%% Cleaning the data by PCA

Y_new = cell{1,n_mat};
En = cell{1,n_mat};
for i = 1: n_mat
    Y = Yn{i};
    Y_bar = Y - mean(Y,2);    
    
    cov_Y = cov(Y_bar');
    [E,evals] = eig(cov_Y);
    [~,idx] = sort(diag(evals),'descend');
    E = E(:,idx(1:Rn));
    
    X = E'*Y_bar;
    Y_new{i} = E*X;
    En{i} = E;
end

% Initializing Zn

Zn = cell{1,n_mat};
P  = zeros(d,m);
D = zeros(d,C);
for i = 1:n_mat
    Zn{i} = randn(Rn);
end

for iter = 1:n_iter
    
    for i = 1:n_mat
        P = P + En{i}*Zn{i};      
    end
    [U,~,V] = svd(P);
    D_new = U(:,1:C)*V(:,1:C)';

    if (sum((D - D_new)^2) <= eta)
    break
    end
    D = D_new;
    for i = 1: n_mat
        Zn{i} = En{i}'*D;
    end
end
    
Xn = cell{1,n_mat};
DXn = cell{1,n_mat};
for i = 1:n_mat
    Xn{i} = (D'*D)\D'*Y_new{i};
    DXn{i} = D*Xn{i};
    DnXn = Y_new{i} - DXn{i};
    
end
    




