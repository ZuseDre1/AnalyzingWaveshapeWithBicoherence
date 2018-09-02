function [cs,coh] = compute_cross_spectrum(data,segleng,segshift,epleng,maxfreqbin)

    % function for the calculation of cross-spectra.
    % 
    % INPUT: 
    % data       : ndat x nchan matrix: each colum is the time-series in one
    %              channel.
    % segleng    : length of each segment in bins, e.g. segleng=1000;  
    % segshift   : numer of bins by which neighboring segments are shifted;
    %              e.g. segshift=segleng/2 makes overlapping segments
    % epleng     : length of each epoch
    % maxfreqbin : max frequency in bins
    %
    % OUTPUT: 
    % cs  : nchan by chan by maxfreqbin by nseg tensor cs(:,:,f,i) contains 
    %       the cross-spectrum at frequency f and segment i
    % coh : complex coherency calculated from cs    
    %
    % USAGE:
    % [cs,coh] = ...
    % compute_cross_spectrum(data,segleng,segshift,epleng,maxfreqbin)
    %
    % This function is a simplified version of a function, which belongs to
    % the METH Toolbox by Guido Nolte:
    %     
    % METH Toolbox: 
    % https://www.uke.de/dateien/institute/neurophysiologie-und-pathophysiologie/downloads/meth.zip
    % (see data2cs_event.m)
    
    zeropad      = 0;
    [ndat,nchan] = size(data);
    maxfreqbin   = min([maxfreqbin,floor((segleng+zeropad)/2)+1]);
    nep          = floor(ndat/epleng);
    nseg         = floor((epleng-segleng)/segshift)+1; 
    cs           = zeros(nchan,nchan,maxfreqbin,nseg); 
    av           = zeros(nchan,maxfreqbin,nseg);
    mywindow     = repmat(hanning(segleng),1,nchan);

    nave = 0;
    % average over epochs;
    for j = 1:nep
        dataep = data((j-1)*epleng+1:j*epleng,:);
        % average over segments;
        for i = 1:nseg 
            dataloc    = dataep((i-1)*segshift+1:(i-1)*segshift+segleng,:);
            datalocfft = fft(dataloc.*mywindow);
            % loop over frequencies
            for f = 1:maxfreqbin     
                cs(:,:,f,i)   = cs(:,:,f,i)+conj(datalocfft(f,:)'*datalocfft(f,:)); 
                av(:,f,i)     = av(:,f,i)+conj(datalocfft(f,:)');                       
            end
        end
        nave = nave+1;
    end
    cs = cs/nave;
    av = av/nave;

    for f = 1:maxfreqbin    
        for i = 1:nseg
            cs(:,:,f,i) = cs(:,:,f,i)-av(:,f,i)*av(:,f,i)';
        end
    end

    ndim = length(size(cs));
    if ndim == 3
        [~,~,n3] = size(cs);
        coh = cs;
        for i = 1:n3
            c          = squeeze(cs(:,:,i));
            coh(:,:,i) = c./sqrt(diag(c)*diag(c)');
        end
    elseif ndim == 4
        [~,~,n3,n4] = size(cs);
        coh           = cs;
        for i = 1:n3
            for j = 1:n4
                c = squeeze(cs(:,:,i,j));
                coh(:,:,i,j) = c./sqrt(diag(c)*diag(c)');
            end
        end
    end
end
    
    
   
    
    
    
    
    
    


