clearvars
close all

% Ef = 20
E = linspace(0,200,21);


for i=1:length(E)
    fprintf('step : %d\n',i)
    ampli(i) = abs(bistability(E(i)));
end

%%
plot(E, ampli)