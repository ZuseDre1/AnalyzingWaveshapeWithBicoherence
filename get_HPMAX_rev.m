function [trans_data,top_filters,top_topos] = get_HPMAX_rev(data,segleng,segshift,epleng,band,maxfreqbins,topMany,dStop,df)

% function for the calculation of a HPMax Spatial Filter
% 
% INPUT: 
% data       : [nrChannels x time] matrix: each row is the time-series in one
%              channel.
% segleng    : length of each segment in bins, e.g. segleng=1000;  
% segshift   : numer of bins by which neighboring segments are shifted;
%              e.g. segshift=segleng/2 makes overlapping segments
% epleng     : length of each epoch
% maxfreqbin : max frequency in bins
% band       : array of frequencies on which HPMax centers the spatial
%              filter
% topMany    : number of spatial filters to be returned by HPMax. Cannot be
%              larger than the number of channels 
% dStop      : frequency bins to estimate the signal power. For a 10 Hz
%              target, df=1 and a frequency resolution of 1 Hz, the signal
%              power is estimated using the power at frequencies 10 +- 1Hz.
%
% df         : number of frequency bins used to estimate the noise power on
%              either site of the target frequency. Depending on the target
%              frequency (and corresponding distance to its harmonics),
%              this parameter should be chosen to ensure that no overlap
%              occurs with harmonics of the target frequency. 
%              For example: for a target frequency f0,
%              df should be chosen such that f0+dStop+df < 2*f0-dStop. 
%              
%              
% OUTPUT: 
% trans_data  : [topMany x time] filterd channel data
% top_filters : [topMany x nrChannels] HPMax spatial filters
% top_topos   : [nrChannels x topMany] Spatial topographies associated with
%               each extracted signal component in trans_data.                
%


Cross_spectrum     = compute_cross_spectrum(data',segleng,segshift,epleng,maxfreqbins);
Cross_spectrum_new = permute(Cross_spectrum,[3,1,2]);
C                  = Cross_spectrum_new;

% extract signal portion for relevant frequencies
C_all                   = real(C(1:maxfreqbins,:,:));
[top_filters,top_topos] = compute_HPMAX_rev(C_all,band,topMany,dStop,df);
trans_data              = cell(size(top_filters));

for i=1:size(top_filters,2)
    cfilt         = top_filters{i};
    trans_data{i} = cfilt*data;
end

