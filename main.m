

%% update Tblock with xls table and EEG data

%%

NrowsXls = size(TblockXls,1); loadedFileName = 'nicnemame'; %loadedFileName = '210330_122311_000';
while true
        % get IDblock if the process has been interrupted
        IDblock = getTableNextRow(Tblock,'done'); % find not processed row
        
%         notDoneID = find(~(Tblock.done(1:IDblock)));
%         if ~isempty(notDoneID)
%             IDblock = min(notDoneID);
%         end
        if IDblock>NrowsXls, break; end;
        
        %file name
        fn = TblockXls.FileName{IDblock};
        
        if ~strcmp(fn,loadedFileName) % if the file we want is different from the one we already have or have anything loaded
            % get path to the file
            fRowInflist = find (strcmp([fn '.smrx'],{flist.name}));  % corresponding row
            fp = [flist(fRowInflist).folder '\' fn '.smrx'] ; % path with file name to the file 

            if length(fRowInflist)>1, disp('Error, nalezeny dva fajly se stejnym nazvem'); return; end;

            % load
            lf = loadSmrxV4(fp);
            fs = lf.T.SamplingFreq(lf.T.Title == 'sig');
            
            downFactor = 2;
            fs = fs/downFactor;
            lf.T = downsampleTable(lf.T,downFactor);
           
            
            
            % file info
            loadedFileName = fn;

            rezNumber = regexpi(fp,'(?<=rez(.*))[123456789]+(?=\\[123456789]+)','match'); % rez, bud mezere nebo ne, pak cislo jakkoliciferny, lomitko a cislo
            if numel(rezNumber)~=1, disp('Error, neidentifikoval jsem cislo rezu'); return; end;
            rezNumber = str2num(rezNumber{1});
            
            inFileSamplePointer = 1; % this will serve us as a pointer that points to the part of the file to be analysed
            
            
%          end
%          
%          if true
            
            %downFactor = 2;
           % lf.T = downsampleTable(lf.T,downFactor);
            % preprocessing
            % filtering
            
            
            
            sig = double(lf.T.Signal{lf.T.Title == 'sig'});
            NsamplesInFile = numel(sig);
            

            s = preprocessInVitroLFP(sig,fs);
            
            clear sig
            
            
            % preprocessing of stimulation signal
            % transform stimulation signal to logicals and compute onsets
            treshold = 2.5;
            stim = double(lf.T.Signal{lf.T.Title == 'opti.stim'});
            stim(1:2)=zeros(1,2); stim(end)=0; % corrections
            stimonofs = getOnsetOffsetByThreshold(stim ,treshold);
            clear stim
            
            s.stim = stimonofs.squareSignal;
            s.stimOnsets = stimonofs.onSignal;
            clear stimonofs
      
            s.t_dn = linspace(lf.timeRecStartDn,lf.timeRecEndDn,NsamplesInFile)'; % timeline in datenum (just for fun...ha ha)
            
            %plot(s.stim)
              
            
            %s.timeRecStartChar = lf.timeRecStartChar;
            
            clear lf
            
            
        end
        
        % we have the correct file loaded, preprocessed, now let's find a block of EEG,
        % analyze it and save.
        Tblock.FileName(IDblock)=fn; % save filename
        Tblock.date(IDblock)=datestr(s.t_dn(1),'dd/mm/yyyy'); % date
        Tblock.slice(IDblock)=num2str(rezNumber); % slice
        Tblock.fs(IDblock)= fs;
        Tblock.IDblock(IDblock)=IDblock; 
        %Tblock.FileName{rp}=
        
        %zpracovat ty spatne bloky
        Tblock.valid(IDblock) = isempty(TblockXls.notes{IDblock});
        
        % zpracovat krizky
        Tblock.light(IDblock) = getMultipleChoice(TblockXls(IDblock,:),labelsLight);
        Tblock.pulse(IDblock) = getMultipleChoice(TblockXls(IDblock,:),labelsPulses);
        disp('analysis');
        
        % napocitat parametry aktualnniho bloku
        %Tblock = getBlockParameters(Tblock);
        
        %sb = getBlockSignal(s)
        
        % get where next stimulation starts 
        stimOnsetInd = find(s.stimOnsets);
        futureStimOnsetsInd = stimOnsetInd(stimOnsetInd>inFileSamplePointer);
        nextStartInd = futureStimOnsetsInd(1);
        nextEndInd = futureStimOnsetsInd(    p.(char(Tblock.pulse(IDblock))).pulses   ) + round(p.(char(Tblock.pulse(IDblock))).periodSec*fs);
        nextLenInd = nextEndInd - nextStartInd +1;
        
        
        
        % check if the observed stimuli matches the claimed one.
        claimedStimDurationSec = p.(char(Tblock.pulse(IDblock))).totalDurationSec;
        claimedDuty = p.(char(Tblock.pulse(IDblock))).percentHi;
        
        deltaLen = (nextEndInd-nextStartInd) - round(claimedStimDurationSec*fs)  ;
        deltaDuty = sum(s.stim(nextStartInd:nextEndInd))/nextLenInd   -  claimedDuty;
        
        if deltaLen>1 | deltaDuty>0.001
            disp('nekonzistence mezi stimula?ním protokolem a skute?n? zm??enými daty!!!!');
        end

      
       Tblock.startSec(IDblock)=(nextStartInd/fs);
       Tblock.endSec(IDblock)=(nextEndInd/fs);
       sc = cropStruct(s,nextStartInd,nextEndInd);
       %plot(Tblock.signal{rp}.noDCnoArt,'k')
       %Tblock.signal{IDblock} = sc;
       %save(fullfile(p_block_signal, [ num2str(IDblock) '.mat']),'sc'); 
       disp('block saved to a table');
       
%% Now we have all recordings in Tblock table, run through all the blocks again but now address the IEDs

       % now lets compute onsets of IEDs from both high frequnency RMS and
       % higest derivative of median (the lower freq component)
       % this part of the script is very cool indeed
        windowN = 200;
       [env,thr] = thresholdTwoClusterSignal(sc.abshf, windowN,false);
       %
       diffMed = [abs(diff(sc.med)); 0];
       windowN = 300;
       [env2,thr2] = thresholdTwoClusterSignal( diffMed , windowN,false);

       windowN = 100;
       [env3,thr3] = thresholdTwoClusterSignal( env2.*env , windowN,false);
       %plot(env3); hold on; plot(thr3*ones(size(env3)));
       
       % set some reasonable margin in ms that every event closer than that
       % will be 
       % connect everything shorter than
       connectMargin_ms = 250;
       connectMargin_N = round(fs*connectMargin_ms/1000);
       envTh = getConsecutiveElements(env3>thr3,connectMargin_N);
       

       of = getOnsetOffsetByThreshold(envTh,thr3);
       sc.IEDonsets = of.onSignal;
       IEDonoffInd = of.onoffdur;
%        plot(sc.noDC);
%        hold on;
%        plot(env3);
%        plot(IEDonsets);
%        plot(0.9*duration);

        % This is so cool, we have onsets of IEDs !!! great work, Emsik.
        % You deserve a medal. And definitely not from president Zeman!

       % now we will put each IED into a IED table

       stimOnsetsInd=find(sc.stimOnsets);
  
       Nied=numel(IEDonoffInd(:,1));
       for kied=1:Nied
     
           %save to table
           IDied = getTableNextRow(Tied,'IDied');
           Tied.IDied(IDied) = IDied; 
           Tied.IDblock(IDied) = IDblock;
           Tied.startInBlockI(IDied) = IEDonoffInd(kied,1);

           allStimsBeforeIED = stimOnsetsInd<IEDonoffInd(kied,1);
           Tied.lastStimBefore(IDied) = IEDonoffInd(kied,1) -  max( stimOnsetsInd(allStimsBeforeIED) );
           % now width , amp, frequency...
           
           % find the end of HFO
           
           preMargin = 120;
           
           Ind=ooRectify(IEDonoffInd(kied,[1 2])+[preMargin 2500],[1 numel(sc.abshf)]);
           
           IEDsignal=sc.noDCnoArt(Ind(1):Ind(2));
           
           
           minN=400;
           maxN = 2000;
           connectShorterThanN = 100;
           
           [oo,stats] = signal2eventsByMinMaxLength(sc.abshf(Ind(1):Ind(2)),minN,maxN,connectShorterThanN);
           oo
           tic
           [oo,stats] = signal2eventsByMinMaxLength(sc.abshf,minN,maxN,connectShorterThanN);
           toc
           
           %%% tohle pujde pryc
           [env,thr] = thresholdTwoClusterSignal(sc.abshf(Ind(1):Ind(2)), 200,false);
           
           iedWidth = min(find(env<thr))+preMargin;
               
           
            Fpass = 200;
            Fstop = 150;
            Apass = 0.5;
            Astop = 65;
         
            d = designfilt('highpassiir', ...
              'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
              'PassbandRipple',Apass,'StopbandAttenuation',Astop, ...
              'DesignMethod','butter','SampleRate',fs);
            %fvtool(d)
            %fpass = 0.05; 
            %[s_hp,dhp] = highpass(lf.T.Signal{lf.T.Title == 'sig'}(a:b),fpass,fs,'ImpulseResponse','iir','Steepness',0.5);  
            
            IEDsignalHFO = filtfilt(d,IEDsignal);
            
            [env,thr] = thresholdTwoClusterSignal(abs(IEDsignalHFO), 1,true);
            
          
           
           Tied.IEDwidth(IDied) = iedWidth/fs;
           
           Ind=ooRectify(IEDonoffInd(kied,[1 2])+[0 500],[1 numel(sc.med)]);
           medMinMax = sc.med(Ind(1):Ind(2));
           
           Tied.IEDamp(IDied) = max(medMinMax)-min(medMinMax);
           
           Ind=ooRectify(IEDonoffInd(kied,[1 2])+[preMargin iedWidth+preMargin],[1 numel(sc.noDCnoArt)]);
           [f,absH,faxis]= getFreqStockwell(sc.noDCnoArt(Ind(1):Ind(2)),fs);
           Tied.HFOfreq(IDied)  = f;
           Tied.HFOpwr(IDied) =   mean(sc.hf(IEDonoffInd(kied,1):IEDonoffInd(kied,2)).^2);   % rms value of filtred signal fd of the size rms_length
           
           plotit;

       end
       
       
       
%        sm = medfilt1(sc.abshf,2000);
%        plot(sm)
%         [pks,locs,w,p] = findpeaks(sc.abshf);
%         
%         rmsLen_ms=3;
%         tic
%         srms2 = getRMSpower(sc.hf,fs,rmsLen_ms);
%         toc
%         
%         % spectrogram approach
%         [sspec,fspc,tspc] = spectrogram(sc.hf,64,62,64,fs);
%         sspec = abs(sspec.^2);
%         
        % vuyziju uz ulozeny signal

        
        disp(['computed ' num2str(IDblock)]);
        Tblock.done(IDblock)=true;
        clear sc
        clear env env2 env3 diffMed
 
end








function y = getMultipleChoice(tb,fields)
% this can be standalone function
% z krizu tabulky udela multiple choice popisek
y = [];
    for i = 1:numel(fields)
        field = fields{i};
        if tb.(field){1}=='x'
            if isempty(y)
                y = field; 
            else
                disp('2 corosses found, faulty input data!!! Erika fix it:D just joking')
                return;
            end
              
        end
    end
    
end


function  y = cropStruct(s,nextStartInd,nextEndInd)
% this can be standalone function
    % crops signal structure
    fn = fieldnames(s);
    for i=1:numel(fn)
        y.(fn{i}) = s.(fn{i})(nextStartInd:nextEndInd); 
    end
end





function s = preprocessInVitroLFP(sig,fs)
% to do : trochu poladit     nastd = 12; nstd = 8;
% this can be stanalone function
            % filter parameters
            Fpass = 0.2;
            Fstop = 0.05;
            Apass = 0.5;
            Astop = 65;
         
            d = designfilt('highpassiir', ...
              'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
              'PassbandRipple',Apass,'StopbandAttenuation',Astop, ...
              'DesignMethod','butter','SampleRate',fs);
            %fvtool(d)
            %fpass = 0.05; 
            %[s_hp,dhp] = highpass(lf.T.Signal{lf.T.Title == 'sig'}(a:b),fpass,fs,'ImpulseResponse','iir','Steepness',0.5);  
            
            noDC = filtfilt(d,sig);

            
            % compute median and high frequency part
            medfiltPoints = round(0.04*fs);
            s.med = medfilt1(noDC,medfiltPoints);
            s.hf =  noDC - s.med;
            
            % now remove artefacts 
            nastd = 12;
            nstd = 8;
            stdNoDC = std(noDC); % std of no dc signal
            medHFlen = round(0.0075*fs); % 50
            medLlen = round(0.0025*fs); %150
            medL = medfilt1(noDC,medLlen);
            s.abshf=abs(s.hf);
            s.medHF = medfilt1(s.abshf,medHFlen);
            thickenFilt =round(fs*5/20000);
            ub = medL + max(s.medHF*nastd,stdNoDC*nstd ); % bounds
            lb = medL - max(s.medHF*nastd,stdNoDC*nstd );
            isOverBounds = noDC>ub | noDC<lb;
            clear ub
            clear lb
            
            s.isOverBounds = filtfilt(ones(thickenFilt ,1),1,double(isOverBounds))>0.5;
            s.noDCnoArt = noDC;
            clear noDC;
            s.noDCnoArt(s.isOverBounds)=medL(s.isOverBounds);
            
            % recompute median and high freqency part for the corrected
            % signal
            %s.med = medfilt1(s.noDC,medfiltPoints);
            s.hf =  s.noDCnoArt - s.med;
            
%             plot(s.noDC);
%             hold on;
%             plot(s.noDCnoArt,'y');
%             plot(s.isOverBounds,'r');
%             plot(ub,'r')
%             plot(lb,'r')

end



function T = downsampleTable(T,fsFactor)
% downsample by fsFactor
%;
for r = 1:size(T,1)
    T.Signal{r} = downsample(T.Signal{r},fsFactor);
    fsOrig = double(T.SamplingFreq(r));
    T.SamplingFreq(r) =  fsOrig / fsFactor;
    
end
            
            
end




function y = getOnsetOffsetByThreshold(x,treshold)
% squareSignal andonSignal offSignal are of length of x
% onoffduration is of length nummber of onsets
% on of duration
% on off duration
            x=x(:);
           y.squareSignal = x > treshold;
           so = diff(y.squareSignal);
           y.onSignal = [so>0;  false ];
           y.offSignal = [so<0;  false ];    
           
           ons = find(y.onSignal);
           offs =  find(y.offSignal);
           y.onoffdur = [ons offs  offs-ons]; 
end


function y = getTableNextRow(T,column)
if iscategorical(T.(column))
    isu = find(isundefined(T.(column)));
    y = isu(1);
end
if isnumeric(T.(column)) || islogical(T.(column))
     isz = find(~(T.(column)));
    y = isz(1);
end
end


