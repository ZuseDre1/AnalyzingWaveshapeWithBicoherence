% Script illustrating relation between biphases and waveshape features 
% of five exemplary periodic signals
% Analyzing the Waveshape of Brain Oscillations with Bicoherence:
% Simulations of
% - Section 1.1
% - Figure 1


close all 
clear
set(0,'defaultfigurecolor',[1 1 1]);

colors   = get(0,'DefaultAxesColorOrder');
allylims = [];
xderlims = [-0.05,0.05;-0.1,0.1;-0.05,0.05;-0.1,0.1;-0.1,0.1];
fs       = 1000;
segleng  = fs;
T        = 60;
mf       = 50;
t        = [0:(1/fs):T-(1/fs)];
turq     = [0,169,157]/255;
black    = [0,0,0];
f0       = 10;
phis     = [0,pi/2,pi,3*pi/2,-3*pi/4]+1.5*pi;

nrHarmonics        = 200;
all_harmonics      = zeros(length(phis),nrHarmonics,length(t));
periodic_signals   = zeros(length(phis),length(t));

skewness1  = zeros(1,5);
dskewness1 = zeros(1,5);
skewness2  = zeros(1,5);
dskewness2 = zeros(1,5);

for phix = 1:length(phis)
    phi0 = phis(phix);
    for h = 1:nrHarmonics
        all_harmonics(phix,h,:)   = (1/(h))*cos(2*pi*h*f0*t+(h-1)*phi0);
    end
    % sum and normalize harmonics
    c_sig                         = sum(squeeze(all_harmonics(phix,:,:)),1);
    ms                            = max(abs(c_sig));
    periodic_signals(phix,:)      = c_sig/ms;
end    

% approximate derivatives of periodic signals
df                    = diff(periodic_signals')';
periodic_signals_drv  = df/(1/fs);
periodic_signals_drv  = [periodic_signals_drv(:,1),periodic_signals_drv];
c_der                 = max(abs(periodic_signals_drv),[],2);
periodic_signals_drv  = periodic_signals_drv./repmat(c_der,1,size(periodic_signals_drv,2));

% compute bicoherences
bicoh_all = cell(1,5);
flat_all  = cell(1,5);

for phix = 1:length(phis)
    mod_signal      = periodic_signals(phix,:);
    mod_signal      = awgn(mod_signal,6);
    [bs,bsnr]       = compute_bispectrum(mod_signal',fs,fs/2,fs,50);
    bicoh           = bs./bsnr;
    bicoh_all{phix} = squeeze(bicoh);
end

%biphases
[X,Y]      = meshgrid([2:mf]);
flat_bicoh = squeeze(bicoh(1,2:mf,2:mf)).';
[xho,yho]  = pol2cart(angle(flat_bicoh(:)),abs(flat_bicoh(:)),4);

% plot configurations
subPlotsX         = 5;   
subPlotsY         = length(phis); 
dataPointsPlot    = 500;
xlimsp1p2         = [0 t(dataPointsPlot)];
ylimsConsistency1 = [-1 1];
ylimsConsistency2 = [-0.2,0.2];
waveSubIdx        = 1;
skew2subIdx       = 2;
asym2subIdx       = 3;
biphase2subIdx    = 5;
pattern2subIdx    = 4;

f = figure;
set(f,'Position',[60 159 984 639]);
set(f,'defaultAxesColorOrder',[black;turq]);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% plot signals and derivatives                                                   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

for phix = 1:length(phis)
    
    subplot(subPlotsY,subPlotsX,(phix-1)*subPlotsX+waveSubIdx)
    
    yyaxis left;
    plot(t(1:dataPointsPlot), periodic_signals(phix,1:dataPointsPlot),'k','LineWidth',1.2)
    hold on
    ylim(ylimsConsistency1);
    skewness1(phix) = skewness(periodic_signals(phix,:),0);
    
    yyaxis right
    phi0=phis(phix);
    plot( t(1:dataPointsPlot), periodic_signals_drv(phix,1:dataPointsPlot),'LineWidth',0.3,'Color',turq)
    dskewness1(phix) = skewness(periodic_signals_drv(phix,:),0);
    ylim(ylimsConsistency2);
    ylim(xderlims(phix,:));
    xlim([0,0.5]) 
    if phix==1
        h = title('Signal','FontWeight','normal','Interpreter','latex','FontSize',12);
        P = get(h,'Position');
    end
    if phix==5
        l=legend('Signal','Derivative');
        pc=get(l,'Position');
        ps=get(gca,'Position');
        set(l,'Position',[ps(1),ps(2)-0.1,pc(3),pc(4)]) 
    end
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% plot signal skewness                                                   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

for phix=1:length(phis)

    subplot(subPlotsY,subPlotsX,(phix-1)*subPlotsX+skew2subIdx)
    phi0 = phis(phix);
    hold off
    [fp,xi]         = ksdensity(periodic_signals(phix,:));
    dummy           = xi;
    skewness2(phix) = skewness(fp);
    
    patch([fp zeros(size(fp))], [dummy flip(dummy)],'r',...
                'FaceColor',colors(phix,:),...
                'LineWidth',0.8,...
                'EdgeColor','none');
    ylim(ylimsConsistency1)
    xlim([0,5])
   
    set(gca,'box','on')
    grid on
    if phix==1
        h = title('Value Distribution','FontWeight','normal','Interpreter','latex','FontSize',12);
        P = get(h,'Position');
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% plot derivative skewness                                                   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
for phix=1:length(phis)
    
    subplot(subPlotsY,subPlotsX,(phix-1)*subPlotsX+asym2subIdx)
    phi0 = phis(phix);
    
    data    = periodic_signals_drv(phix,:);
    [fp,xi] = ksdensity(data);
    dummy   = xi;
    patch([fp zeros(size(fp))], [dummy flip(dummy)],'r',...
                   'FaceColor', colors(phix,:),...
                   'LineWidth', 0.8,...
                   'EdgeColor', 'none');
    dskewness2(phix) = skewness(fp);
    ylim(ylimsConsistency2)
    ylim(xderlims(phix,:));
	set(gca,'box','on')
    grid on
    xlim([0,50])
    if phix==1
        h = title('Derivative Distribution','FontWeight','normal','Interpreter','latex','FontSize',12);
        P = get(h,'Position');
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% plot bicoherence patterns                                                 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
for phix=1:length(phis)
    subplot(subPlotsY,subPlotsX,(phix-1)*subPlotsX+pattern2subIdx)
    bb = abs(bicoh_all{phix});
    imagesc(bb)
    caxis([0,1])
    %cptcmap('viridis','mapping','scaled','flip',false);
    
    if phix==1
        h = title('Bicoherence Pattern','FontWeight','normal','Interpreter','latex','FontSize',12);
        P = get(h,'Position');        
    end
    
    if phix==5
       c=colorbar('southoutside');
       pc=get(c,'Position');
       ps=get(gca,'Position');
       set(c,'Position',[ps(1),ps(2)-0.15,pc(3)-0.5*pc(3),pc(4)+0.5*pc(4)]) 
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% plot biphases                                                  
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
for phix=1:length(phis)
    subplot(subPlotsY,subPlotsX,(phix-1)*subPlotsX+biphase2subIdx)
    gp       = get(gca,'Position');
    gp(1)    = gp(1)-0.07*gp(1);  
    biphases = (0.9*bicoh_all{phix});
    c        = compass(biphases(f0+1,f0+1));
    
    set(c,'Color',colors(phix,:))
    set(c,'LineWidth',1.5)
    set(gca,'Position',gp)
    if phix==1
        h = title('Biphase','FontWeight','normal','Interpreter','latex','FontSize',12);
        P = get(h,'Position');
    end
end
set(findall(gcf, 'String', '30', '-or','String','60', '-or','String','150', '-or','String','120','-or','String','0.5') ,'String', '  ');
set(findall(gcf, 'String', '210', '-or','String','240', '-or','String', '300', '-or','String','330') ,'String', '  ');
set(findall(gcf, 'String', '0'),'String', '0','FontSize',12);
set(findall(gcf, 'String', '90'),'String', '\pi/2','FontSize',12);
set(findall(gcf, 'String', '180'),'String', '\pi','FontSize',12);
set(findall(gcf, 'String', '270'),'String', ' 3\pi/2','FontSize',12);
set(findall(gcf, 'String', '  0.5', '-or','String','1'),'String', ''); 
set(findall(gcf, 'String', '  1', '-or','String','1'),'String', ''); 




