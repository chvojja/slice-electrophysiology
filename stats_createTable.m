
% plot options

% alpha = 0.65;
% velikostKolecka = 80;
tloustkaCary = 1.2;


%% generate data
Tblock = Tblock(1:76,:);
Tied = Tied(1:4181,:);
% 
% 
 T = innerjoin(Tblock,Tied,'LeftKeys',1,'RightKeys',2);

 T.G = findgroups(T.FileName,T.slice); % slice grouping

%%
% filterID=T.light=='noLight' | T.light=='light__' | T.pulse =='x500ms' | (T.pulse =='x100ms' &  );
% Tf=T(filterID,:);
% 
% Tf.G_date_slice = findgroups(Tf.FileName,Tf.slice);

Ns=numel(unique(T.G));


% 
% vv = varfun(@mean,T,'InputVariables','IEDamp',...
%        'GroupingVariables','G');

st = struct;
plotDims = [17 22];

Ntests = 4;
% 
% Ttes = table(['M';'F';'M'],[45 45;41 32;40 34],...
%     {'NY';'CA';'MA'},[true;false;false]);

% writetable(T,fXLSname,'Sheet','width','Range','D2');


%%
featureNames = {'HFOwidth_ms','HFOfreq','HFOpwr','IEDamp'};
units = {'ms','Hz','mv^2','mv'};
pulseType = 'x100ms';
%treatmentLight = {'light__','light___'};
treatmentLight = {'light___'};


fXLSname = ['D:\Google Drive - st47058\#PhD#Analysis\2021-10-22 git erika slices\' pulseType '.xlsx'];

for i = 1:Ntests
    featureName = featureNames{i};
    st = getStats(T, st,Ns, featureName,pulseType,treatmentLight);
    
    out = st.(featureName);
    d = out(:,2)-out(:,1);
    percents = 100*([out(:,2)./out(:,1)]-1);
    [d_mean, d_SEM] = meansem(d,3);
    [percents_mean, percents_SEM] = meansem(percents,3);
    [ctrl_mean,ctrl_SEM] = meansem(out(:,1),3);
    [exp_mean,exp_SEM] = meansem(out(:,2),3);
    
    
%     st.TprismPercent(i,:) = [zeros(1,11) percents'  ];
%     %st.TprismRatios(i,:) = [ones(1,11) [out(:,2)./out(:,1)]'  ];
%     st.Tprism(i,:) = [out(:,1)'  out(:,2)'];
    [p(i),h] = signrank(percents);
    % [h, p(i)] = ttest(percents);
%     if i ==3
%         [p(i),h] = signrank(percents)
%     else
%     [h, p(i)] = ttest(percents);
%     end
    
    % Levy sloupec
    writecell({['Paired observations of ' featureName ' in units: ' units{i}]},fXLSname,'Sheet',featureName,'Range','G1'); % header
    writecell({'control','treatment','difference','percent change'},fXLSname,'Sheet',featureName,'Range','C2'); % header
    writematrix([out, d, percents],fXLSname,'Sheet',featureName,'Range','C3');
    %writetable(T,fXLSname,'Sheet',featureName,'Range','B3');
    writematrix([1:numel(d)]',fXLSname,'Sheet',featureName,'Range','B3'); % header
    
    writecell({'mean';'SEM'},fXLSname,'Sheet',featureName,'Range','A15'); % header
    writematrix([ctrl_mean exp_mean],fXLSname,'Sheet',featureName,'Range','C15'); % header
    writematrix([ctrl_SEM exp_SEM],fXLSname,'Sheet',featureName,'Range','C16'); % header
    
    
    % Pravy soupec
    writecell({'p value, uncorrected';'p value, Holm-Sidak correction'},fXLSname,'Sheet',featureName,'Range','J2'); % header 2
    writematrix(p(i),fXLSname,'Sheet',featureName,'Range','K2');
    
    % mean sem mean sem
    writecell({'mean of perc.change';'SEM of perc.change'},fXLSname,'Sheet',featureName,'Range','J5'); % header 
    writematrix([percents_mean; percents_SEM],fXLSname,'Sheet',featureName,'Range','K5');
    writecell({'mean of diff';'SEM of diff'},fXLSname,'Sheet',featureName,'Range','J7'); % header 
    writematrix([d_mean; d_SEM],fXLSname,'Sheet',featureName,'Range','K7');
end


[c_pvalues] = fwer_sidak(p, 0.05);

for i = 1:Ntests
    featureName = featureNames{i};
    writematrix(c_pvalues(i),fXLSname,'Sheet',featureName,'Range','K3');
end






%% functions


function st = getStats(T, st,Ns ,featureName,pulseType,treatmentLight)

st.(featureName)=NaN(Ns,2); % first row is control, second is  tratment
st.groups = [];
dPoint=1;
for groupNum=1:Ns

Ts = T(T.G==groupNum,:);
% Ts.light,Ts.pulse
[G_light_pulse,GTsID] = findgroups(Ts(:,[8 9]));
%G_light_pulse = findgroups(Ts.light,Ts.pulse);
Ts.G = G_light_pulse;

controlID = find(GTsID.light == 'noLight' & GTsID.pulse ==pulseType);
treatmentID = find(GTsID.light == treatmentLight & GTsID.pulse ==pulseType);

if (~isempty(controlID)) & (~isempty(treatmentID) )
Y = varfun(@mean,Ts,'InputVariables',featureName,'GroupingVariables','G');
st.(featureName)(dPoint,1) = Y{controlID,3};
st.(featureName)(dPoint,2) = Y{treatmentID,3};
st.groups(end+1)=groupNum;
dPoint=dPoint+1;  
end


end
st.(featureName)(dPoint:end,:)=[]; % delete empty
%st.(featureName)(:,1);

[st.p,st.h] = signrank(st.(featureName)(:,1),st.(featureName)(:,2));
%st.controlMean=st.(featureName)(:,1);
end




function [means,sems] = getAllHistCounts(dataset,subjects,edges)
% dataset has each column for each subject and may contain NaNs
perctls = [0.01 99.99];
%perctls = [0.05 99.5];
    for col = subjects
    notNaNs=~isnan(dataset(:,col));
    data = dataset(notNaNs,col);
    % tds=st.timeDiffLvsCms(notNaNs,col);
    % tds = rmoutliers(tds,'percentiles',perctls);
    % histogram(tds,edges,'Normalization','probability');
    % hold on;
    % tds=st.timeDiffLvsPLms(notNaNs,col);
    % tds = rmoutliers(tds,'percentiles',perctls);
    % histogram(tds,edges,'Normalization','probability');


    x = rmoutliers(data,'percentiles',perctls);

    [hCounts,edges,bin] = histcounts(x,edges,'Normalization','probability');

    hCountsAll(col,:)=hCounts;

    end
 N=numel(subjects);
 if N>1
    means= mean(hCountsAll(subjects,:));
    sems=std(hCountsAll(subjects,:))/sqrt(N);
 else
    means =  hCountsAll(subjects,:);
    sems=zeros(size(means));
 end
end


function plotMeanSEM(means,sems,velikostKolecka,tloustkaCary,color,alpha,typ)
% means=[10 5 8 4 3];
% sems=[2 1 1.5 0.8 0.5 ];
% timeline = 0.5:1:4.5;
% velikostKolecka = 150;
% tloustkaCary = 2;

Nbars=numel(means);

for i=1:Nbars
    x=i;
    y = means(i);
    err = sems(i);
    
    switch nargin
        case 7
        hb = scatter(x,y,velikostKolecka,'ko','LineWidth',tloustkaCary); 
        set(hb,'MarkerFaceColor',color);  
        
        
        case 6
        hb = bar(x,y,'LineWidth',tloustkaCary,'FaceColor',color,'FaceAlpha',alpha);
   
    end
        
    
    hold on;
    errorbar(x,y,-err,err,'LineWidth',tloustkaCary,'Color',color);
    
end




set(gca,'XLim',[0 Nbars+1],'YLim',[0 max(means+sems)*1.1]);
xticks =1:Nbars;
set(gca,'XTick',xticks);
set(gca,'TickLength',[0 0]);
% XTickLabel {'a','b'...}
hold off;
end


function plotBoxPlot(data,velikostKolecka,alpha,filledon)
color = [0 0 0];
wspacing=2;
dabs = 0.4;
Ncol = size(data,2);
xbars = [1:Ncol]*wspacing-dabs;
xpoints = [1:Ncol]*wspacing+dabs;

if size(data,1)>100
boxplot(data, 'positions', xbars, 'labels', xbars,'Colors','k','PlotStyle','traditional','symbol', '','Notch','on');
else
boxplot(data, 'positions', xbars, 'labels', xbars,'Colors','k','PlotStyle','traditional','symbol', '');
end;

hold on;

for col=1:Ncol
    scdata = data(~isnan(data(:,col)),col);
    rozmitani = ([(randn(size(scdata))-0.5)/25]);
    hb = scatter(xpoints(col)*ones(size(scdata))+rozmitani,scdata,velikostKolecka); %,'MarkerEdgeColor',color
    
    if nargin >3
        hb = scatter(xpoints(col)*ones(size(scdata))+rozmitani,scdata,velikostKolecka,'filled'); 
    end
    set(hb,'MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerEdgeAlpha',alpha,'MarkerFaceAlpha',alpha);
    
end

plot(xpoints(:,[1 2]),data')
set(gca,'XLim',[0 max(xpoints) + 2*dabs]);
box off;
end

