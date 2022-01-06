% plotting the IED 


%% when nothing is plotted
if ~exist('pt') % init
    % matrix
    pt.Rows=6;
    pt.Cols=6;
    pt.TotSubs = pt.Rows*pt.Cols;
    pt.subPlotPointer = 0;
    figure('color','white');
    
end


% if (pt.subPlotPointer == 1) % if first subplot to be plotted
%     figure;
% 
% end



%% Plotting each IED
pt.subPlotPointer=pt.subPlotPointer+1;
subplot(pt.Rows,pt.Cols,pt.subPlotPointer);

IEDsig = sc.noDCnoArt(IEDInd(1):IEDInd(2));
N_IEDsig=numel(IEDsig);
tIED = linspace(0,1000*N_IEDsig/fs,N_IEDsig);

plot(tIED,IEDsig,'k'); hold on;

xLow = abs(IEDInd(1)-HFOInd(1));
xHigh = (xLow + HFOwidthInd);

ylowHigh = ooStretchPercent([min(IEDsig) max(IEDsig)],[1.2 1.2]);
ooPlotPatches(1000*[xLow xHigh]/fs,ylowHigh,'r',0.5);

box off;
%axis on;
ax=gca;
ax.YAxis.Visible = 'off';
xlim([0 numel(IEDsig)]*1000/fs);
xlabel('time, ms');



% 
% subplot(pt.Rows,pt.Cols,pt.subPlotPointer);
% oo = IEDonoffInd(kied,[1 2]);
% IndIED=ooRectify(oo+[preMargin HFOwidthInd+preMargin],[1 numel(sc.noDCnoArt)]);
% IndIEDwider=ooRectify(oo+[0 preMargin + 2500],[1 numel(sc.noDCnoArt)]);
% widerSig = sc.noDCnoArt(IndIEDwider(1):IndIEDwider(2));
% plot(widerSig,'k'); hold on;
% 
% xLow = abs(IndIED(1)-IndIEDwider(1));
% xHigh = xLow + abs(IndIED(1)-IndIED(2));
% 
% ylowHigh = ooStretchPercent([min(widerSig) max(widerSig)],[1.2 1.2]);
% ooPlotPatches([xLow xHigh],ylowHigh,'r',0.6);
% 
% box off;
% axis off;



% Plot Spectrogram
pt.subPlotPointer=pt.subPlotPointer+1;
subplot(pt.Rows,pt.Cols,pt.subPlotPointer);

imagesc(absH);
NfreqTicsk=6;
Nfre = numel(faxis);
yticklabels = round(linspace(faxis(1),faxis(end),NfreqTicsk));
yticks = linspace(1, Nfre, numel(yticklabels));
set(gca, 'YTick', round(yticks), 'YTickLabel', flip(yticklabels(:)));
set(gca,'XTick',[])
ylabel('frequency, Hz');
xlabel(['IED_ID:' num2str(IDied) ' IDblock:' num2str(IDblock)]);

hold on;
redLineY = Nfre-round(Nfre*(f-faxis(1))/abs(faxis(end)-faxis(1)));
line(1:size(absH,2),redLineY*ones(1,size(absH,2)),'LineWidth',1.2,'Color','r');

% HFO detail
pt.subPlotPointer=pt.subPlotPointer+1;
subplot(pt.Rows,pt.Cols,pt.subPlotPointer);
hfowInd=round((HFOInd(2)-HFOInd(1))/4);
HFOdetail = sc.noDCnoArt(HFOInd(1):HFOInd(2));
N_HFOsig=numel(HFOdetail);
tHFO= linspace(0,1000*N_HFOsig/fs,N_HFOsig);
plot(tHFO,HFOdetail,'k'); hold on;

box off;
%axis on;
ax=gca;
ax.YAxis.Visible = 'off';
xlim([0 numel(HFOdetail)]*1000/fs);
xlabel('time, ms');

if pt.subPlotPointer == pt.TotSubs % figure is full;
     hf = gcf;
    fName = [ 'D:\tempErikaSST\plots_IEDs'   '\lastIDblock_'   num2str(IDblock) '_lastIDied_' num2str(IDied)  ];
   %  fName = [ 'D:\tempErikaSST\plots_IEDs'   '\lastIDblock_'  ];
   
    print2png(hf,fName,[30  22],80);
    pause(0.100);

     close(gcf);
    clear pt;
end