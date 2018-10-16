function [ yi ] = interpmean( x, y, xi )
%INTERPMEAN interpolation (downsampling) by averaging
%   YI = INTERPMEAN(X, Y, XI) returns the values of the 1D function (X, Y)
%   between points XI (with typically length(XI)<<length(X)) by taking the
%   nanmean of Y on each interval [XI(i) XI(i+1)]

Ni = length(xi);
yi = nan(1, Ni-1);
for i_sample = 1:Ni-1
    %fprintf(['interpol ' num2str(i_sample) ' of ' num2str(Ni-1) '\n']);
    idx = x>=xi(i_sample) & x<=xi(i_sample+1);
    yi(i_sample)=nanmean(y(idx));
end

end
