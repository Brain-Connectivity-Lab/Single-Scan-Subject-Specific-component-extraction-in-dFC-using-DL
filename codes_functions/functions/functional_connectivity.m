function vec = functional_connectivity(reg_ts,Z,win,stride)

% reg_ts Timeseries of each ROI   T x ROI x sub
% win  window size for DFC
% stride for DFC
[T,ROI,sub] = size(reg_ts);

if nargin < 4
    stride = 1;
if nargin < 3
    win = T;
    stride = 1;
    if nargin < 2
       Z = 0;
    end
end
end
cor = zeros(ROI);
if Z == 1
    msg1 = 'With Fisher Z transform ..';
else 
    msg1 = '..';
end
if win == T
    msg = ['Creating Static map ',msg1];
else
    msg = ['Creating Dynamic map ',msg1];
end
x = 0;
f = waitbar(x,msg);
for s = 1:sub
    p=1;    
   for i1 = 1:stride:T-win+1
        cor(:,:,s) = corr(reg_ts(i1:i1+win-1,:,s));

    if nnz(isnan(cor(:,:,s)))
        disp('nan in')
        disp(s)
    end
    if Z == 1
    cor(:,:,s) = 0.5*(log(1+cor(:,:,s))-log(1-cor(:,:,s)));
    end
    vec(:,p,s) = corrvec(cor(:,:,s));
    p=p+1;
   end
    x = s/size(reg_ts,3);
    waitbar(x,f) 
end
close(f)
vec = squeeze(vec);
end