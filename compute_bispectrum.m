                          
function [bs,bsnr] = compute_bispectrum(data,segleng,segshift,epleng,maxfreqbins)

% function for the calculation of bispectra and bicoherence.  
%
% INPUT: 
% data        : ndat times nchan matrix each colum is the time-series in one
%               channel;
% segleng     : length of each segment in bins, e.g. segleng=1000;  
% segshift    : numer of bins by which neighboring segments are shifted;
%               e.g. segshift=segleng/2 makes overlapping segments23333
% epleng      : length of each epoch
% maxfreqbins : maximum frequency in bins, starting at zeros Hertz. 
%               The frequency resolution, df, is given by the physical length of a
%               segment, say T. Then df=1/T. E.g. if T=2 seconds, the maxfreqbins=101
%               means that the maximum physical frequency is 50 Hertz.
%
% OUTPUT: 
% bs   : [nchan x nf x nf] tensor for nf frequencies (i.e. nf = maxfreqbins)   
%        bs(i,f1,f2)=<x(f1)_i*x(f2)_i*conj(x(f1+f2-1)_i)>
%        where  x is the Fourier-transform of the data of each segment
%
% bsnr : corresponding normalization factor defined by 
%        bsn(i,f1,f2)=N_i(f1) N_i(f2) N_i(f1+f2-1);
%        where N_p(f) is defined as (<abs(x(f)_p)^3>)^(1/3) 
%        Bicoherence can be calculated as bs./bsn 
%
% nave : number of averages
%
% USAGE: 
% [bs,bsnr]=compute_bispectrum(data,segleng,segshift,epleng,maxfreqbins);
%
% AUTHOR: 
% This function is a simplified version of a function, which belongs to
% the METH Toolbox by Guido Nolte:
%     
% METH Toolbox: 
% https://www.uke.de/dateien/institute/neurophysiologie-und-pathophysiologie/downloads/meth.zip
% (see data2bs_univar.m)



[ndat,nchan] = size(data);
nf           = maxfreqbins;
mywindow     = repmat(hanning(segleng),1,nchan);

nep  = floor(ndat/epleng);  
nseg = floor((epleng-segleng)/segshift)+1;
bs   = zeros(nchan,nf,nf);
bsnr = zeros(nchan,nf,nf);
bsn  = zeros(nchan,2*nf-1);
nave = 0;

% for each episode
for j = 1:nep   
    dataep = data((j-1)*epleng+1:j*epleng,:); 
    % average over all segments
    for i = 1:nseg 
        dataloc    = dataep((i-1)*segshift+1:(i-1)*segshift+segleng,:);
        datalocfft = fft(detrend(dataloc).*mywindow); %compute windowed fft of cut out segment
        datalocfft = datalocfft(1:2*nf-1,:);   
        bslocn     = ((abs(datalocfft)).^3)';
        % for each channel
        for ichan = 1:nchan 
          xx               = hankel(conj(datalocfft(1:2*nf-1,ichan)));
          bsloc(ichan,:,:) = (datalocfft(1:nf,ichan)*transpose(datalocfft(1:nf,ichan))).*xx(1:nf,1:nf);
          bs(ichan,:,:)    = bs(ichan,:,:)+bsloc(ichan,:,:);
        end
        nave = nave+1;
        bsn  = bsn+bslocn;
    end
end



bs  = bs/nave;
bsn = bsn/nave;
bsn = power(bsn,1/3);

for i = 1:nchan
    for f1 = 1:nf
        for f2 = 1:nf
            bsnr(i,f1,f2) = (bsn(i,f1)*bsn(i,f2)*bsn(i,f1+f2-1));
        end
    end
end



    

    
