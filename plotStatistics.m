
% % plot options
% 
% alpha = 0.65;
% velikostKolecka = 80;
% tloustkaCary = 1.2;
% 
% 
% %% generate data
% Tblock = Tblock(1:76,:);
% Tied = Tied(1:4181,:);
% 
% 
% T = innerjoin(Tblock,Tied,'LeftKeys',1,'RightKeys',2);
% 
% T.G = findgroups(T.FileName,T.slice);

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

%%
featureName = 'HFOwidth_ms';
pulseType = 'x500ms';
%pulseType = 'x100ms';
treatmentLight = 'light__';

st = getStats(T, st,Ns, featureName,pulseType,treatmentLight);

%%
alpha = 1; velikostKolecka = 12;
plotBoxPlot(st.(featureName),velikostKolecka,alpha);
if st.h
    title(['Differnece in ' featureName ' under pulse ' pulseType ' and ' treatmentLight ' was significant p=' num2str(st.p,2) ]);
    %title(['Differnece in ' featureName ' under pulse ']);
else
    title(['Differnece in ' featureName ' under pulse ' pulseType ' and ' treatmentLight ' was not significant p=' num2str(st.p,2) ]);
end
%set(gca,'XTickLabel',chNames(sh));
hf = gcf;
print2png(hf,[featureName '_' pulseType '_' treatmentLight],plotDims,80);
close(hf);
 
%%
featureName = 'HFOwidth_ms';
pulseType = 'x500ms';
%pulseType = 'x100ms';
treatmentLight = 'light___';

st = getStats(T, st,Ns, featureName,pulseType,treatmentLight);

%%
alpha = 1; velikostKolecka = 12;
plotBoxPlot(st.(featureName),velikostKolecka,alpha);
if st.h
    title(['Differnece in ' featureName ' under pulse ' pulseType ' and ' treatmentLight ' was significant p=' num2str(st.p,2) ]);
    %title(['Differnece in ' featureName ' under pulse ']);
else
    title(['Differnece in ' featureName ' under pulse ' pulseType ' and ' treatmentLight ' was not significant p=' num2str(st.p,2) ]);
end
%set(gca,'XTickLabel',chNames(sh));
hf = gcf;
print2png(hf,[featureName '_' pulseType '_' treatmentLight],plotDims,80);
close(hf);
 






%%
featureName = 'HFOfreq';
pulseType = 'x500ms';
%pulseType = 'x100ms';
treatmentLight = 'light__';

st = getStats(T, st,Ns, featureName,pulseType,treatmentLight);

%%
alpha = 1; velikostKolecka = 12;
plotBoxPlot(st.(featureName),velikostKolecka,alpha);
if st.h
    title(['Differnece in ' featureName ' under pulse ' pulseType ' and ' treatmentLight ' was significant p=' num2str(st.p,2) ]);
    %title(['Differnece in ' featureName ' under pulse ']);
else
    title(['Differnece in ' featureName ' under pulse ' pulseType ' and ' treatmentLight ' was not significant p=' num2str(st.p,2) ]);
end
%set(gca,'XTickLabel',chNames(sh));
hf = gcf;
print2png(hf,[featureName '_' pulseType '_' treatmentLight],plotDims,80);
close(hf);


%%
featureName = 'HFOfreq';
pulseType = 'x500ms';
%pulseType = 'x100ms';
treatmentLight = 'light___';

st = getStats(T, st,Ns, featureName,pulseType,treatmentLight);

%%
alpha = 1; velikostKolecka = 12;
plotBoxPlot(st.(featureName),velikostKolecka,alpha);
if st.h
    title(['Differnece in ' featureName ' under pulse ' pulseType ' and ' treatmentLight ' was significant p=' num2str(st.p,2) ]);
    %title(['Differnece in ' featureName ' under pulse ']);
else
    title(['Differnece in ' featureName ' under pulse ' pulseType ' and ' treatmentLight ' was not significant p=' num2str(st.p,2) ]);
end
%set(gca,'XTickLabel',chNames(sh));
hf = gcf;
print2png(hf,[featureName '_' pulseType '_' treatmentLight],plotDims,80);
close(hf);







%%
featureName = 'HFOpwr';
pulseType = 'x500ms';
%pulseType = 'x100ms';
treatmentLight = 'light__';

st = getStats(T, st,Ns, featureName,pulseType,treatmentLight);

%%
alpha = 1; velikostKolecka = 12;
plotBoxPlot(st.(featureName),velikostKolecka,alpha);
if st.h
    title(['Differnece in ' featureName ' under pulse ' pulseType ' and ' treatmentLight ' was significant p=' num2str(st.p,2) ]);
    %title(['Differnece in ' featureName ' under pulse ']);
else
    title(['Differnece in ' featureName ' under pulse ' pulseType ' and ' treatmentLight ' was not significant p=' num2str(st.p,2) ]);
end
%set(gca,'XTickLabel',chNames(sh));
hf = gcf;
print2png(hf,[featureName '_' pulseType '_' treatmentLight],plotDims,80);
close(hf);




%%
featureName = 'HFOpwr';
pulseType = 'x500ms';
%pulseType = 'x100ms';
treatmentLight = 'light___';

st = getStats(T, st,Ns, featureName,pulseType,treatmentLight);

%%
alpha = 1; velikostKolecka = 12;
plotBoxPlot(st.(featureName),velikostKolecka,alpha);
if st.h
    title(['Differnece in ' featureName ' under pulse ' pulseType ' and ' treatmentLight ' was significant p=' num2str(st.p,2) ]);
    %title(['Differnece in ' featureName ' under pulse ']);
else
    title(['Differnece in ' featureName ' under pulse ' pulseType ' and ' treatmentLight ' was not significant p=' num2str(st.p,2) ]);
end
%set(gca,'XTickLabel',chNames(sh));
hf = gcf;
print2png(hf,[featureName '_' pulseType '_' treatmentLight],plotDims,80);
close(hf);





%%
featureName = 'IEDamp';
pulseType = 'x500ms';
%pulseType = 'x100ms';
treatmentLight = 'light__';

st = getStats(T, st,Ns, featureName,pulseType,treatmentLight);

%%
alpha = 1; velikostKolecka = 12;
plotBoxPlot(st.(featureName),velikostKolecka,alpha);
if st.h
    title(['Differnece in ' featureName ' under pulse ' pulseType ' and ' treatmentLight ' was significant p=' num2str(st.p,2) ]);
    %title(['Differnece in ' featureName ' under pulse ']);
else
    title(['Differnece in ' featureName ' under pulse ' pulseType ' and ' treatmentLight ' was not significant p=' num2str(st.p,2) ]);
end
%set(gca,'XTickLabel',chNames(sh));
hf = gcf;
print2png(hf,[featureName '_' pulseType '_' treatmentLight],plotDims,80);
close(hf);



%%
featureName = 'IEDamp';
pulseType = 'x500ms';
%pulseType = 'x100ms';
treatmentLight = 'light___';

st = getStats(T, st,Ns, featureName,pulseType,treatmentLight);

%%
alpha = 1; velikostKolecka = 12;
plotBoxPlot(st.(featureName),velikostKolecka,alpha);
if st.h
    title(['Differnece in ' featureName ' under pulse ' pulseType ' and ' treatmentLight ' was significant p=' num2str(st.p,2) ]);
    %title(['Differnece in ' featureName ' under pulse ']);
else
    title(['Differnece in ' featureName ' under pulse ' pulseType ' and ' treatmentLight ' was not significant p=' num2str(st.p,2) ]);
end
%set(gca,'XTickLabel',chNames(sh));
hf = gcf;
print2png(hf,[featureName '_' pulseType '_' treatmentLight],plotDims,80);
close(hf);





%% functions

function st = getStats(T, st,Ns ,featureName,pulseType,treatmentLight)

st.(featureName)=NaN(Ns,2); % first row is control, second is  tratment

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
st.(featureName)(groupNum,1) = Y{controlID,3};
st.(featureName)(groupNum,2) = Y{treatmentID,3};
end

end

%st.(featureName)(:,1);

[st.p,st.h] = signrank(st.(featureName)(:,1),st.(featureName)(:,2));
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

