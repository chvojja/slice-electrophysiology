% plotting the IED 


%% when nothing is plotted
if ~exist('pt') % init
    % matrix
    pt.Rows=6;
    pt.Cols=8;
    pt.TotSubs = pt.Rows*pt.Cols;
    pt.subPlotPointer = 1;
    figure('color','white');
    
end


% if (pt.subPlotPointer == 1) % if first subplot to be plotted
%     figure;
% 
% end



%% Plotting each IED

subplot(pt.Rows,pt.Cols,pt.subPlotPointer);
oo = IEDonoffInd(kied,[1 2]);
IndIED=ooRectify(oo+[preMargin iedWidth+preMargin],[1 numel(sc.noDCnoArt)]);
IndIEDwider=ooRectify(oo+[0 preMargin + 2500],[1 numel(sc.noDCnoArt)]);
widerSig = sc.noDCnoArt(IndIEDwider(1):IndIEDwider(2));
plot(widerSig,'k'); hold on;

xLow = abs(IndIED(1)-IndIEDwider(1));
xHigh = xLow + abs(IndIED(1)-IndIED(2));

ylowHigh = ooStretchPercent([min(widerSig) max(widerSig)],[1.2 1.2]);
ooPlotPatches([xLow xHigh],ylowHigh,'r',0.6);

box off;
axis off


pt.subPlotPointer=pt.subPlotPointer+1;



if pt.subPlotPointer == pt.TotSubs % figure is full;
     hf = gcf;
    fName = [ 'D:\tempErikaSST\plots_IEDs'   '\lastIDblock_'   num2str(IDblock) '_lastIDied_' num2str(IDied)  ];
   %  fName = [ 'D:\tempErikaSST\plots_IEDs'   '\lastIDblock_'  ];
   
    print2png(hf,fName,[30  22],80);
   % ;
%     wh = [18  16];
%     quality = 90;
%     
%     set(hf,'renderer','painters');
%     set(hf, 'Units', 'normalized');
%     set(hf, 'Position',[0 0 wh(1) wh(2)]);
% 
%     set(hf, 'PaperPositionMode', 'manual', ...
%             'PaperPosition', [0 0 wh(1) wh(2)], ...
%             'InvertHardCopy', 'on');
%     print( 'sdsds.jpg'); 
%     
%     
%     close(gcf);
    
end