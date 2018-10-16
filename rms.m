function [ out ] = rms( in )
%RMS root mean squared

out = sqrt(mean(in.^2));

end

