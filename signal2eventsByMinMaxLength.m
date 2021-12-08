function [oo,stats] = signal2eventsByMinMaxLength(sig,minN,maxN,connectShorterThanN)
%SIGNAL2EVENTSBYMINMAXLENGTH Summary of this function goes here
%   Detailed explanation goes here
%   clusters timeserie signal based on whatever into two distinct clusters
% minN maxN minimum and maximum size in samples of the events.

numberOfVariables = 1;
xSorted=sort(sig,'Ascend');
lb = [xSorted(1)];
ub = [xSorted(end)];

% ranges = [lb ub];

fun = @(x) fitness(x,sig,minN,maxN,connectShorterThanN);

% [thrOut,fval] = ga(fun,numberOfVariables,[],[],[],[],lb,ub); % genetic je
% pomalej

%x0 = lb + rand(size(lb)).*(ub - lb);
 x0=lb+0.5*(ub - lb);
 %x0=0.015
 
 options = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf);
 
 %[thrOut,fval] = patternsearch(fun,x0,[],[],[],[],lb,ub,[], options) ; % tohle je cool rychly
% [thrOut,fval] = patternsearch(fun,x0,[],[],[],[],lb,ub) ; % tohle je cool rychly
[thrOut,fval] = surrogateopt(fun,lb,ub)


squareSignal = thresholdSignal(thrOut,sig,minN,maxN,connectShorterThanN);



so = diff(squareSignal);
onSignal = [so>0;  false ];
offSignal = [so<0;  false ];    

ons = find(onSignal);
offs =  find(offSignal);
oo = [ons offs  offs-ons]; 

stats.thr=thrOut;
stats.squareSignal=squareSignal;


    function y = fitness(x,sig,minN,maxN,connectShorterThanN)
    % fitness function 
    % x bude prah ktery optimalizuju
    
    % get threshold
    tsig = thresholdSignal(x,sig,minN,maxN,connectShorterThanN);

%      plot(sig); hold on;
%     plot(x*ones(size(sig)),'r');
    
    % plot(tsig,'b');
    
    % assess it
    
    if isempty(tsig) % pokud to nevybere nic
        y = Inf;
    else
        
        ups=mean(sig(tsig))+var(sig(tsig));
        lows=mean(sig(~tsig))+var(sig(~tsig));

        y = 1/((ups-lows).^2); %y =1-ofInBounds
    end

    end


    function tsig = thresholdSignal(x,sig,minN,maxN,connectShorterThanN)
        
        sig=sig(:); % sig is column
        tsig = sig>x;
        tsig(1)=false; tsig(end) = false;

  

        J = imclose(tsig,ones(connectShorterThanN,1)); % spojovani, zapln mezery mensi nez structural element
        J = imopen(J,ones(minN,1)); % odstran sum/eventy kratsi  jak structural element. vsechno co je mensi nez [1 1 1] odstrani a nacha vsechny [1 1 1] a delší atd takze odstrani sum mensi nez 3 puntiky
        if numel(find(J))<2
            tsig=[];
            return
        end
        invJ = imclose(~J,ones(maxN,1)); % odstran vetsi nebo rovno nez
        tsig = J & invJ; 
    end


end
