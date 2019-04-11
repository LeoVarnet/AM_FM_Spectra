function [ FM_withnan ] = remove_artifacts_FM( FM_withnan, fs, maxjump, minduration, frange, varargin )
%REMOVE_ARTIFACTS_FM Remove the artifacts from f0 or FM trajectories. 
% Use as [ FM_no_artifacts ] = remove_artifact_FM( FM, fs, maxjump, minduration, frange )
%   maxjump: maximum gap in frequency between two samples
%   minduration: durée minimum d'un segment isolé
%   frange: limits of credible pitch range (frange = [minf maxf])
%
% Three steps: - suppress jumps between consecutive samples larger than maxjump
%              - supress segments shorter than minduration
%              - supress extreme values (< minf or > maxf)
%
% The new version of the function (old version : 'remove_artifact_FM')
% includes several additional modules, corresponding to optional
% parameters:
% [ FM_no_artifacts ] = remove_artifact_FM( __ , frange_median, var_thres, plot_step ) 
%    frange_median: limits of credible pitch range, relative to the median (= [minf_median maxf_median])
%    var_thres: removal of fricative with a threshold of (average) variability of FM
%    plot_step: before/after plot
% 
% Optional steps: - supress extreme values relative to the median (< minf_median*median or > maxf_median*median)
%   (deprecated:) - supress values outside of the "global range" of the sound

% Leo Varnet - 2017 (last modified 04/2019)

minf = frange(1); maxf = frange(2);
if ~isempty(varargin)
    frange_median = varargin{1};
end
if length(varargin)>1
    var_thres = varargin{2};
else
    var_thres = 0;
end
if length(varargin)>2
    plot_step = varargin{3};
else
    plot_step = 'no';
end

if isyes(plot_step)
    figure;
    plot((1:length(FM_withnan))/fs, FM_withnan, 'r'); xlabel('time (s)'); ylabel('frequency (Hz)'); ylim([0 800]); hold on
end

% suppress jumps between consecutive samples larger than maxjump
deriveFM=diff(FM_withnan);
FM_withnan(abs(deriveFM)>maxjump)=NaN;

% supress segments shorter than minduration
isFM = ~isnan(FM_withnan);
i=1;
nminsamples=minduration*fs;

while i<length(isFM)
    if isFM(i)==1
        startpoint=i;
        while i<length(isFM) & isFM(i)==1
            i=i+1;
        end
        endpoint=i;
        if endpoint-startpoint<nminsamples
            FM_withnan(startpoint:endpoint)=NaN;
        elseif var_thres>0
        % calculate the average variability
            segment_var=nansum(abs(diff(FM_withnan(startpoint:endpoint))))/((endpoint-startpoint)/fs);
            if segment_var>var_thres
                FM_withnan(startpoint:endpoint)=NaN;
            end
        end
    end
    i=i+1;
end

% supress outliers
% - Absolute values (frange)
FM_withnan(FM_withnan<minf | FM_withnan>maxf ) = NaN;
if isyes(plot_step)
    xlimits = xlim; xlim(xlimits);
    plot(xlimits, [minf minf], 'y--'); plot(xlimits, [maxf maxf], 'y--');
end
% - Relative values (frange_median)
if ~isempty(varargin)    
    medianFM = nanmedian(FM_withnan);
    FM_withnan(FM_withnan<frange_median(1)*medianFM | FM_withnan>frange_median(2)*medianFM ) = NaN;
    if isyes(plot_step)
        xlimits = xlim; xlim(xlimits);
        plot(xlimits, [1 1]*frange_median(1)*medianFM, 'c--'); plot(xlimits, [1 1]*frange_median(2)*medianFM, 'c--'); plot(xlimits, [1 1]*medianFM, 'c:');
    end
end

% [n,x]=hist(FM_withnan,50);
% n = n(2:end-1); x = x(2:end-1);
% [~,idxFMmode]=max(n);
% i=idxFMmode;
% while i<length(n) & n(i)>0
%     i=i+1;
% end
% maxfbis=x(i)+(x(2)-x(1));
% i=idxFMmode;
% while i>1 & n(i)>0
%     i=i-1;
% end
% minfbis=x(i)-(x(2)-x(1));
% FM_withnan(FM_withnan<minfbis | FM_withnan>maxfbis )=NaN;
% 

if isyes(plot_step)
    plot((1:length(FM_withnan))/fs, FM_withnan, 'g');
end

end

