function [oo,stats] = signal2eventsByMinMaxLength(sig,minN,maxN,connectCloserThanN)
%SIGNAL2EVENTSBYMINMAXLENGTH Summary of this function goes here
%   Detailed explanation goes here
%   clusters timeserie signal based on whatever into two distinct clusters
% minN maxN minimum and maximum size in samples of the events.

xSorted=sort(sig,'Ascend');
lb = [xSorted(1)];
ub = [xSorted(end)];

% fun = @(x) fitness(x,sig,maxN);
%  
errorY=0.001;
errorX=10e-8; % tahle presnost je klicova

% check if lower bound adjustment is required:
if numel(sig)>maxN  
    if 1==adjustLowerBound(lb,sig,maxN,connectCloserThanN) % if we break the maxN rule even at lb
     
        lb = bisectionMethod( @(x)adjustLowerBound(x,sig,maxN,connectCloserThanN) ,lb,ub,errorY,errorX);
    end
end
% adjust upper bound

ub = bisectionMethod( @(x)adjustUpperBound(x,sig,minN,connectCloserThanN) ,lb,ub,errorY,errorX);

% now we will optimize inside lb and ub

% lb
% ub

% lb = 0.0055;
% ub = 0.044;
% 
% lb = 0.003;
% ub = 0.05;

fun = @(x) fitness(x,sig,minN,connectCloserThanN);
%options = optimset('OutputFcn', @outfun);
% [x, fval] = fminbnd(fun, lb, ub, options);
[thr, fval] = fminbnd(fun, lb, ub);

%  options = optimoptions('surrogateopt','MaxFunctionEvaluations',400);
%  [thr, fval] = surrogateopt(fun,lb,ub,options);


tout = thresholdSignal(thr,sig,minN,connectCloserThanN);



so = diff(tout);
onSignal = [so>0;  false ];
offSignal = [so<0;  false ];    

ons = find(onSignal);
offs =  find(offSignal);
oo = [ons offs  offs-ons]; 

stats.thr=thr;
stats.tsig=tout;


% plot(sig); hold on;
% plot(thr*ones(size(sig)),'r');
% plot(tout,'b');
% fval


%disp('done');
%% optimalizacni funkce
    function y = fitness(x,sig,minN,connectCloserThanN)
    % fitness function 
    % x bude prah ktery optimalizuju
    
    % get threshold
    tsig = thresholdSignal(x,sig,minN,connectCloserThanN);


        ups=mean(sig(tsig)); %+var(sig(tsig));
        lows=mean(sig(~tsig)); %+var(sig(~tsig));
        
%         ups=var(sig(tsig));
%         lows=var(sig(~tsig));

        %y = 1/((ups-lows).^2);
        y = min(ups,lows)/max(ups,lows);
        if isnan(ups) ||  isnan(lows)
            y = Inf;
        end


    end


    function tsig = thresholdSignal(x,sig,minN,connectCloserThanN)
        
        sig=sig(:); % sig is column
        tsig = sig>x;
        tsig(1)=false; tsig(end) = false;


        tsig = imclose(tsig,ones(connectCloserThanN,1)); % spojovani, zapln mezery mensi nez structural element
        tsig = imopen(tsig,ones(minN,1)); % odstran sum/eventy kratsi  jak structural element. vsechno co je mensi nez [1 1 1] odstrani a nacha vsechny [1 1 1] a delší atd takze odstrani sum mensi nez 3 puntiky

% %         
%         plot(sig); hold on;
%         plot(x*ones(size(sig)),'r');
%         plot(tsig,'b');
    end


%% podpurne funkce pro lb ub bounds
    function y = adjustLowerBound(x,sig,maxN,connectCloserThanN)
        sig=sig(:); % sig is column
        tsig = sig>x;
        tsig(1)=false; tsig(end) = false;
        
        tsig = imclose(tsig,ones(connectCloserThanN,1)); % spojovani, zapln mezery mensi nez structural element
        
       [~,len] = logi2ooInd(tsig);
       NoutBounds = numel(find(len>maxN));  % kriteriem je zanik vetsiho eventu nez je maxN
       if NoutBounds>0
           y =1;
       else
           y =-1;
       end

    end 

    function y = adjustUpperBound(x,sig,minN,connectCloserThanN)
        sig=sig(:); % sig is column
        tsig = sig>x;
        tsig(1)=false; tsig(end) = false;
        
        tsig = imclose(tsig,ones(connectCloserThanN,1)); % spojovani, zapln mezery mensi nez structural element
        
       [~,len] = logi2ooInd(tsig);
       NlongerThanMinN = numel(find(len>minN)); % kriterium je  po spojeni mensich vznik alespon jedne o minimalni velikosti minN
       if NlongerThanMinN>0
           y =-1;
       else
           y =1;
       end

    end 

%     function stop = outfun(x,optimValues,state)
%         stop = false;
%         % Check whether directional derivative norm is less than .01.
%         if norm(optimValues.directionalderivative) < .01
%             stop = true;
%         end 
%     end

end

% 
% function [oo,stats] = signal2eventsByMinMaxLength(sig,minN,maxN,connectShorterThanN)
% %SIGNAL2EVENTSBYMINMAXLENGTH Summary of this function goes here
% %   Detailed explanation goes here
% %   clusters timeserie signal based on whatever into two distinct clusters
% % minN maxN minimum and maximum size in samples of the events.
% 
% numberOfVariables = 1;
% xSorted=sort(sig,'Ascend');
% lb = [xSorted(1)];
% ub = [xSorted(end)];
% 
% % ranges = [lb ub];
% 
% fun = @(x) fitness(x,sig,minN,maxN,connectShorterThanN);
% 
% % [thrOut,fval] = ga(fun,numberOfVariables,[],[],[],[],lb,ub); % genetic je
% % pomalej
% 
% %x0 = lb + rand(size(lb)).*(ub - lb);
%  x0=lb+0.5*(ub - lb);
%  %x0=0.015
%  
%  options = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf);
%  
%  %[thrOut,fval] = patternsearch(fun,x0,[],[],[],[],lb,ub,[], options) ; % tohle je cool rychly
% % [thrOut,fval] = patternsearch(fun,x0,[],[],[],[],lb,ub) ; % tohle je cool rychly
% [thrOut,fval] = surrogateopt(fun,lb,ub)
% 
% 
% squareSignal = thresholdSignal(thrOut,sig,minN,maxN,connectShorterThanN);
% 
% 
% 
% so = diff(squareSignal);
% onSignal = [so>0;  false ];
% offSignal = [so<0;  false ];    
% 
% ons = find(onSignal);
% offs =  find(offSignal);
% oo = [ons offs  offs-ons]; 
% 
% stats.thr=thrOut;
% stats.squareSignal=squareSignal;
% 
% 
%     function y = fitness(x,sig,minN,maxN,connectShorterThanN)
%     % fitness function 
%     % x bude prah ktery optimalizuju
%     
%     % get threshold
%     tsig = thresholdSignal(x,sig,minN,maxN,connectShorterThanN);
% 
% %      plot(sig); hold on;
% %     plot(x*ones(size(sig)),'r');
%     
%     % plot(tsig,'b');
%     
%     % assess it
%     
%     if isempty(tsig) % pokud to nevybere nic
%         y = Inf;
%     else
%         
%         ups=mean(sig(tsig))+var(sig(tsig));
%         lows=mean(sig(~tsig))+var(sig(~tsig));
% 
%         y = 1/((ups-lows).^2); % + 1-ofInBounds;
%     end
% 
%     end
% 
% 
%     function tsig = thresholdSignal(x,sig,minN,maxN,connectShorterThanN)
%         
%         sig=sig(:); % sig is column
%         tsig = sig>x;
%         tsig(1)=false; tsig(end) = false;
% 
%   
% 
%         J = imclose(tsig,ones(connectShorterThanN,1)); % spojovani, zapln mezery mensi nez structural element
%         J = imopen(J,ones(minN,1)); % odstran sum/eventy kratsi  jak structural element. vsechno co je mensi nez [1 1 1] odstrani a nacha vsechny [1 1 1] a delší atd takze odstrani sum mensi nez 3 puntiky
%         if numel(find(J))<2
%             tsig=[];
%             return
%         end
%        % invJ = imclose(~J,ones(maxN,1)); % odstran vetsi nebo rovno nez
%        % tsig = J & invJ; 
%     end
% 
% 
% end
