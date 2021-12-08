function [oo,stats] = signal2eventsByMinMaxLength(sig,minN,maxN,connectShorterThanN)
%SIGNAL2EVENTSBYMINMAXLENGTH Summary of this function goes here
%   Detailed explanation goes here
%   clusters timeserie signal based on whatever into two distinct clusters
% minN maxN minimum and maximum size in samples of the events.

FitnessFunction = @(x) fitness(x,sig,minN,maxN,connectShorterThanN);
numberOfVariables = 1;
xSorted=sort(sig,'Ascend');

lb = [xSorted(1),xSorted(1)];
ub = [xSorted(end),xSorted(end)];
[out,fval] = ga(FitnessFunction,numberOfVariables,[],[],[],[],lb,ub);


disp('');
oo=out;
stats = [];


function y = fitness(x,sig,minN,maxN,connectShorterThanN)
% fitness function for GA
% x bude prah ktery optimalizuju
tsig = sig>x;
tsig(1)=false; tsig(end) = false;

J = imclose(tsig,ones(1,connectShorterThanN)); % spojovani, zapln mezery mensi nez structural element
J = imopen(J,ones(1,minN)); % odstran sum/eventy kratsi  jak structural element. vsechno co je mensi nez [1 1 1] odstrani a nacha vsechny [1 1 1] a delší atd takze odstrani sum mensi nez 3 puntiky
invJ = imclose(~J,ones(1,maxN)); % odstran vetsi nebo rovno nez
tsig = J & invJ;


ups=median(xsig(tsig))+var(xsig(tsig));
lows=median(xsig(~tsig))+var(xsig(~tsig));

y = 1/((ups-lows).^2);
  
end

end
