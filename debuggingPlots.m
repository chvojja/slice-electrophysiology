plot(sc.noDCnoArt)
hold on;


stem(sc.stimOnsets,'r');
stem(sc.IEDonsets,'k');

plot(50*env3,'g');

hold on;
plot(0.01*sc.noDCnoArt)

plot(0.01*envTh,'r')