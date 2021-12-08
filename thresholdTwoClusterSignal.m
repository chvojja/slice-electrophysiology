function [envHigh,thr] =  thresholdTwoClusterSignal(x, windowN, verbose)
% estimates best threshold based on some noisy rectified (power) data
% windowN = delka s jakou se to ma vyhladit
% this is fucking hard function I dont know how to describe it
       [envHigh, envLow] = envelope(x,windowN,'rms');
       [hc,edg]=histcounts(envHigh);
       
       %hc=20*log10(hc);
       %hc(hc==-Inf)=0;  % this is very important@!!!!!!!!
       hc =  medfilt1(hc,ceil(numel(hc)/20));
       hc = [ 0  hc 0];
       %figure; plot( )
       
       [~,locs,~,pp] = findpeaks(hc);
       amplifiedPP=hc(locs).*pp;
       
       [~,sI] = sort(amplifiedPP,'descend');
       firstSecondLocIin_hc = locs(sI(1:2));
       amplifiedPPs=amplifiedPP(sI);
       
       %[sIs,~] = sort(sI(1:2),'ascend');
       [firstSecondLocIin_hcs,~] = sort(firstSecondLocIin_hc,'ascend');
       
       [~,locs2,~,pp2] = findpeaks(-hc(firstSecondLocIin_hcs(1):firstSecondLocIin_hcs(2)));
       amplifiedPPInv=hc(locs2).*pp2;
       
       [sm,smI] = sort(amplifiedPPInv,'descend');
       %midI = locs2(smI(1));
       %locmidI = locs(sIs(1) + smI(1));
       locminimumIin_hc = firstSecondLocIin_hcs(1)+locs2(smI(1))-1;
%       
  
        
% %     
       
     %  if abs(firstSecondLocIin_hc+1)>sum(abs(firstSecondLocIin_hc-locminimumIin_hc))
           thr = edg(locminimumIin_hc);
           
%        else
%            thr = mean(edg(firstSecondLocIin_hc));
%        end

if verbose
        subplot(2,1,1);
        hold on;
        plot(hc);
        plot(locminimumIin_hc,hc(locminimumIin_hc),'ko');
        plot(firstSecondLocIin_hc,hc(firstSecondLocIin_hc),'kx');
        
        subplot(2,1,2);
        plot(envHigh); hold on;
        plot(thr*ones(size(envHigh)),'r');
end
end