function [ AM, FM ] = hilbert_extraction( D, fs, varargin )
% [ AM, FM ] = HILBERT_EXTRACTION( D, fs, thres )
%   AM and FM extraction of D (sampling frequency fs) with Hilbert
%   transform. thres is an optional argument for artifact rejection
%   (consider FM only in regions where AM>thres). 
%
% Leo Varnet 2016 - last modified 11/10/2018

nbch = size(D,1);
if length(varargin)==0
    thres = 0;
elseif length(varargin)==1
    thres = varargin{1};
elseif length(varargin)>1
    error('too many arguments in function hilbert_extraction')
end

t=(1:length(D))/fs;

if isreal(D)
    hilbert_responses = hilbert(D);
else
    hilbert_responses = D;
end

% AM extraction
for k=1:nbch,
    AM(k,:)=abs(hilbert_responses(k,:));
end

if nargout>1
    % FM extraction
    for k=1:nbch
        TFS(k,:) = unwrap(angle(hilbert_responses(k,:)));
        FM(k,:) = (1/(2*pi))*gradient(TFS(k,:),t);
    end
end

for k=1:nbch,
    if thres>0
        FM(k,AM(k,:)<thres)=NaN;
    end
end

end

