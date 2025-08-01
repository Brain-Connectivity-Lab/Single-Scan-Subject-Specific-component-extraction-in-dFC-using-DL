


function [ D, Xn, DnXn, DXn ] = cobe_DL( Yn ,C,nor) 

%% COBE Dictionary Learning Algorithm

% formulation -->        Yn = DXn + DnXn'
% input -->    Yn = a cell containing all the matrices --> size(Yn) = 1 x n (n = number of matrices) 
%              C = no. of atoms in common dictionary D
%              if C == 0 then the function uses cobe algorithm. 
%              nor = 1      zscore normalization on Yn 
%                  = 0      no normalization (Default)

% output -->   D = Common dictionary --> size = m x C (m is the dimention)
if nargin == 2
    nor = 0;
end

A = length (Yn);
for i = 1 : length (Yn)
    Yn{i} = Yn{i}';
end
%%if Z Normalisation is required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if nor == 1
    Yn_bar = cell (A, 1);

    for i = 1: length (Yn)
        X =Yn {1,i};
        A = zscore (X);
        A = A';
        Yn_bar{i,1} = A;
    end
    Yn = Yn_bar;
else
    Yn_bar = cell (A, 1);

    for i = 1: length (Yn)
        X =Yn {1,i};
        A = X;
        Yn_bar{i,1} = A;
    end
    Yn = Yn_bar;
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Getting the components to be removed and the applying COBE
rC =  C;
opts.c = rC;
if C ~=0
[D, Xn] = cobec(Yn, opts);
else 
[D, Xn] = cobe(Yn);
end    
% getting the common values
DXn= cell(1,1);
for i= 1: length(Xn)
    DXn{i,1} = D* Xn{1,i};
end

DnXn= cell(1,1);
for i= 1: length (DXn)
    DnXn{i,1} = Yn{i,1}-DXn{i,1};
end

end

