function [ FM_withnan ] = remove_artifact_FM( FM_withnan, fs, maxjump, minduration, minf, maxf )
%REMOVE_ARTIFACT_FM Remove the artifacts from f0 of FM trajectories. 
% Use as [ f0_no_artifacts ] = remove_artifact_FM( f0, fs, maxjump, minduration, minf, maxf )
%   maxjump : maximum gap in frequency between two samples
%   minduration : durée minimum d'un segment isolé
%   minf and maxf : limits of credible pitch range
% Four steps : - suppress jumps between consecutive samples larger than maxjump
%              - supress segments shorter than minduration
%              - supress extreme values (< minf or > maxf)
%              - supress values outside of the "global range" of the sound

% suppress jumps between consecutive samples larger than maxjump
derivef0=diff(FM_withnan);
FM_withnan(abs(derivef0)>maxjump)=NaN;

% supress segments shorter than minduration
isf0 = ~isnan(FM_withnan);
i=1;
nminsamples=minduration*fs;

while i<length(isf0)
    if isf0(i)==1
        startpoint=i;
        while i<length(isf0) & isf0(i)==1
            i=i+1;
        end
        endpoint=i;
        if endpoint-startpoint<nminsamples
            FM_withnan(startpoint:endpoint)=NaN;
        end
    end
    i=i+1;
end

% supress implausible values
FM_withnan(FM_withnan<minf | FM_withnan>maxf )=NaN;

% [n,x]=hist(FM_withnan,50);
% n = n(2:end-1); x = x(2:end-1);
% [~,idxf0mode]=max(n);
% i=idxf0mode;
% while i<length(n) & n(i)>0
%     i=i+1;
% end
% maxfbis=x(i)+(x(2)-x(1));
% i=idxf0mode;
% while i>1 & n(i)>0
%     i=i-1;
% end
% minfbis=x(i)-(x(2)-x(1));
% FM_withnan(FM_withnan<minfbis | FM_withnan>maxfbis )=NaN;
% 

end

