clearvars
close all

% Ef = 20
dE = 0.5;
Ei = 0;
Ef = 20;
E = Ei:dE:Ef;


for i=1:length(E)
    fprintf('step : %d\n',i)
    ampli(i) = abs(bistability(E(i)));
end

%%
scatter(abs(E), ampli)