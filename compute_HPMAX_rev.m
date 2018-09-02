function [top_evecs,top_topos]=compute_HPMAX_rev(C_all,band,topMany,dStop,df)
    

    nr_f      = size(band,2);
    top_evecs = cell(1,nr_f);
    top_topos = cell(1,nr_f);

    for i = 1:nr_f
    
        apix = band(i);
        wd   = df;           
        C    = 0;
        C    = C+(1/(2*dStop+1))*sum(C_all([apix-dStop:apix+dStop],:,:),1);
        C    = C+(1/(2*dStop+1))*sum(C_all([apix*2-1-dStop:apix*2-1+dStop],:,:),1);
        C    = C+(1/(2*dStop+1))*sum(C_all([apix*3-2-dStop:apix*3-2+dStop],:,:),1);
        C    = squeeze(C);

        D11  = (1/(wd))*sum(C_all([apix-dStop-wd:apix-dStop-1],:,:),1);
        D12  = (1/(wd))*sum(C_all([apix+dStop+1:apix+dStop+wd],:,:),1);
        D1   = 0.5*D11+0.5*D12;
        D21  = (1/(wd))*sum(C_all([apix*2-1-dStop-wd:apix*2-1-dStop-1],:,:),1);
        D22  = (1/(wd))*sum(C_all([apix*2-1+dStop+1:apix*2-1+dStop+wd],:,:),1);
        D2   = 0.5*D21+0.5*D22;
        D31  = (1/(wd))*sum(C_all([apix*3-2-dStop-wd:apix*3-2-dStop-1],:,:),1);
        D32  = (1/(wd))*sum(C_all([apix*3-2+dStop+1:apix*3-2+dStop+wd],:,:),1);
        D3   = 0.5*D31+0.5*D32;
        D    = (D1+D2+D3);


        
        D            = squeeze(D);
        DH           = D^(0.5);
        M            = inv(DH)*C*inv(DH);
        [tvec,tval]  = eig(M);
        [~,I]        = sort(diag(tval),'descend');
        tvec         = tvec(:,I);
        vecs         = inv(DH)*tvec;        
        Evecs        = normc(vecs);        
        top_evecs{i} = Evecs(:,1:topMany)';
        top_topos{i} = D*Evecs(:,1:topMany);
    end
end

