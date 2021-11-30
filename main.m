% %clear all;
% 
% %% Definition of available stimuli
% % labels of stimulation parameters
% labelsLight = {'noLight','light_','light__','light___'};
% labelsPulses = {'x10ms','x500ms','x1s','x3s','x10s'};
% 
% % stimulation protocol structure
% p = struct; 
% p = setProtocolStructByNameWidthPeriodPulses(p,'x10ms',0.01,3,50); % POZOR JE TO 100ms
% p = setProtocolStructByNameWidthPeriodPulses(p,'x500ms',0.5,5,100); %500ms 5s period 100pulses
% p = setProtocolStructByNameWidthPeriodPulses(p,'x1s',1,5,100);
% p = setProtocolStructByNameWidthPeriodPulses(p,'x3s',1,5,5);
% p = setProtocolStructByNameWidthPeriodPulses(p,'x10s',10,17.8,1);
% 
% 
% 
% %% get Erikas .xlsx with information about each block of stimulation 
% fp_Tstim='Z:\Petranova\SST, CHR2 H134R\table_slicestim_SST_CHR2_H134R.xlsx';
% 
% % ImportOptions_ = detectImportOptions(fp_Tstim, 'NumHeaderLines', 1);
% % TblockXls = readtable(fp_Tstim,ImportOptions_);
% % save('TblockXls','TblockXls'); % save for testing purpose
% 
% load('TblockXls','TblockXls');
% 
% % get xls table statistics
% 
% Tstats = table('Size',[5 4],'VariableTypes',repmat({'double'},1,4),'VariableNames',labelsLight,'RowNames',labelsPulses);
% 
% for ip=1:numel(labelsPulses)
%     pulseType = labelsPulses{ip};
%     for il=1:numel(labelsLight)
%         lightType = labelsLight{il};
%         contains = strcmpi(TblockXls.(pulseType),'x') & strcmpi(TblockXls.(lightType),'x') & strcmpi(TblockXls.notes,''); % only "x" and only with empty "note" column
%         Tstats{pulseType,lightType}  = numel(find(contains));
%     end
%     
% end
% 
% 
% %isempty(TblockXls.light_{1})
% 
% 
% %% get a list with all .smrx files and path to them
% p_EEG = 'Z:\Petranova\SST, CHR2 H134R\recordings';
% flist = dir([p_EEG '\**\*.smrx']);  % file list
% 
% 
% 
% Tblock = getTblockTemplate();
% %% update Tblock with xls table and EEG data



NrowsXls = size(TblockXls,1); loadedFileName = 'nicnemame'; %loadedFileName = '210330_122311_000';
for IDblock =24:NrowsXls
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
            fs = fs/2;
            downFactor = 2;
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
            
            s.t_dn = linspace(lf.timeRecStartDn,lf.timeRecEndDn,NsamplesInFile)'; % timeline in datenum (just for fun...ha ha)
            
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
      
            
            %plot(s.stim)
              
            
            %s.timeRecStartChar = lf.timeRecStartChar;
            
            clear lf
            
            
        end
        
        % we have the correct file loaded, preprocessed, now let's find a block of EEG,
        % analyze it and save.
        Tblock.FileName(IDblock)=fn; % save filename
        %Tblock.date(rp)=; % date
        Tblock.slice(IDblock)=num2str(rezNumber); % slice
        Tblock.IDblock(IDblock)=num2str(IDblock); % slice % tohle nevim jestli by nebylo lepsi krokovat od jedne pro kazdy rez
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

      
       sc = cropStruct(s,nextStartInd,nextEndInd);
       %plot(Tblock.signal{rp}.noDCnoArt,'k')
       Tblock.signal{IDblock} = sc;
       disp('block saved to a table');
       
      
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
 
end


%% Now we have all recordings in Tblock table, run through all the blocks again but now address the IEDs
Tied = getTiedTemplate(5000);
for IDblock =1:NrowsXls

       % now lets compute onsets of IEDs from both high frequnency RMS and
       % higest derivative of median (the lower freq component)
       % this part of the script is very cool indeed
        windowN = 200;
       [env,thr] = thresholdTwoClusterSignal(sc.abshf, windowN);
       %
       diffMed = [abs(diff(sc.med)); 0];
       windowN = 300;
       [env2,thr2] = thresholdTwoClusterSignal( diffMed , windowN);
   
       windowN = 100;
       [env3,thr3] = thresholdTwoClusterSignal( env2.*env , windowN);
       plot(env3); hold on; plot(thr3*ones(size(env3)));
       
       
       of = getOnsetOffsetByThreshold(env3,thr3);
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
           IEDonInd(kied);
            
           %save to table
           IDied = getTableNextEmptyRow(Tied,'IDied');
           Tied.IDied(IDied) = IDied; 
           Tied.IDblock(IDied) = IDblock;
           Tied.startInBlockI(IDied) = IEDonoffInd(kied,1);
           
           allStimsBeforeIED = stimOnsetsInd<IEDonoffInd(kied,1);
           Tied.lastStimBefore(IDied) = IEDonoffInd(kied,1) -  max( stimOnsetsInd(allStimsBeforeIED) );
           
           
           
       end
end



%% Support functions
function Tblock = getTblockTemplate()
% this can be standalone function
% definition of a table of blocks
% this table is the main data type in the analysis
% a row in the table represents a "block" - a part of EEG with parameters:
varTypes = {'categorical','categorical','categorical','categorical','cell',  'double' ,'logical','categorical', 'categorical',    'double',  'double'  };
varNames = {'IDblock', 'FileName',     'date',       'slice'    , 'signal','fs',     'valid',       'light',       'pulse',      'startDn','endDn'   };
Tblock = table('Size',[1,numel(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
end

function Tied = getTiedTemplate(Nrows)
% this can be standalone function
% definition of a table of blocks
% this table is the main data type in the analysis
% a row in the table represents a "block" - a part of EEG with parameters:
varTypes = {'categorical', 'categorical', 'double',      'double',      'double',           'double',  'double', 'double'  , 'double'};
varNames = {'IDied',        'IDblock',    'startInBlockI', 'endInBlockI', 'lastStimBefore',   'SLEwidth','SLEamp',  'HFOfreq','HFOpwr'};
Tied = table('Size',[Nrows,numel(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
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


function p_structure = setProtocolStructByNameWidthPeriodPulses(p_structure, name,widthSec,periodSec,pulses)
        % define stimulation sequence
        p_structure.(name).pulses = pulses;
        p_structure.(name).widthSec = widthSec;
        p_structure.(name).periodSec = periodSec;
        p_structure.(name).totalDurationSec = pulses*periodSec; % total length of stimualtion 
        p_structure.(name).percentHi=(pulses*widthSec)/(pulses*periodSec); % percent of time being hi state
         
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
            medHFlen = 0.0075*fs; % 50
            medLlen = 0.0025*fs; %150
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


function rms = getRMSpower(sig,fs,rmsLen_ms)
        rmsLen=round(rmsLen_ms*10^-3*fs);
        rms=zeros(size(sig));

        seg_down=floor((rmsLen-1)/2);
        seg_up=ceil((rmsLen-1)/2);
        for i=1:size(sig,1)
            seg_start=i-seg_down; %nastaveni zacatku segmentu
            seg_end=i+seg_up; %nastaveni konce segmentu
            if seg_start<1 %korekce u zacatku signalu
                seg_start=1;
            end
            if seg_end>size(sig,1) %korekce u konce signalu
                seg_end=size(sig,1);
            end
            rms(i,:)=sqrt(  mean(sig(seg_start:seg_end,:).^2));   % rms value of filtred signal fd of the size rms_length 
        end
        
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

function [envHigh,thr] =  thresholdTwoClusterSignal(x, windowN)
% estimates best threshold based on some noisy rectified (power) data
% windowN = delka s jakou se to ma vyhladit
% this is fucking hard function I dont know how to describe it
       [envHigh, envLow] = envelope(x,windowN,'rms');
       [hc,edg]=histcounts(envHigh);
       
       hc=20*log10(hc);
       hc(hc==-Inf)=0;  % this is very important@!!!!!!!!
       hc =  medfilt1(hc,ceil(numel(hc)/20));
       
       %figure; plot( )
       
       [~,locs,~,pp] = findpeaks(hc);
       amplifiedPP=hc(locs).*pp;
       
       [~,sI] = sort(amplifiedPP,'descend');
       firstSecondLocIin_hc = locs(sI(1:2));
       
       
       %[sIs,~] = sort(sI(1:2),'ascend');
       [firstSecondLocIin_hcs,~] = sort(firstSecondLocIin_hc,'ascend');
       
       [~,locs2,~,pp2] = findpeaks(-hc(firstSecondLocIin_hcs(1):firstSecondLocIin_hcs(2)));
       amplifiedPPInv=hc(locs2).*pp2;
       
       [sm,smI] = sort(amplifiedPPInv,'descend');
       %midI = locs2(smI(1));
       %locmidI = locs(sIs(1) + smI(1));
       locminimumIin_hc = firstSecondLocIin_hcs(1)+locs2(smI(1))-1;
       
%         hold on;
%         plot(hc);
%         plot(locminimumIin_hc,hc(locminimumIin_hc),'ko');
%         plot(firstSecondLocIin_hc,hc(firstSecondLocIin_hc),'kx');
% %     
%        
       if abs(firstSecondLocIin_hc+1)>sum(abs(firstSecondLocIin_hc-locminimumIin_hc))
           thr = edg(locminimumIin_hc);
           
       else
           thr = mean(edg(firstSecondLocIin_hc));
       end
end


function y = getOnsetOffsetByThreshold(x,treshold)
% squareSignal andonSignal offSignal are of length of x
% onoffduration is of length nummber of onsets
% on of duration
% on off duration
           y.squareSignal = x > treshold;
           so = diff(y.squareSignal);
           y.onSignal = [so>0;  false ];
           y.offSignal = [so<0;  false ];    
           
           ons = find(y.onSignal);
           offs =  find(y.offSignal);
           y.onoffdur = [ons' offs'  offs'-ons']; 
end


function y = getTableNextEmptyRow(T,column)
isu = isundefined(T.(column));
y = isu(1);
end