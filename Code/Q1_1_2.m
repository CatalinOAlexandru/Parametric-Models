% NOTE:
% Run Q1_1 before running this file.

Avox = dwis(:,92,65,72);

% Define a starting point for the non-linear fit
startx = [3.5e+00 3e-03 2.5e-01 0 0];

% Define various options for the non-linear fitting algorithm

h=optimset('MaxFunEvals',20000,... 
    'Algorithm','quasi-newton',... 
    'TolX',1e-10,... 
    'Display','off',...
    'TolFun',1e-10);

% Now run the fitting
[parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD',startx,h,Avox,bvals,qhat);

disp(RESNORM);
disp(parameter_hat);
% plot(Avox,RESNORM, 'o');

[dif, Ypred] = BallStickSSD(parameter_hat,Avox,bvals,qhat);
plotting(Ypred, Avox)



function graph = plotting(Ypred, Avox)
    graph = figure();
    plot(Avox, 'rx')
    hold on;
    plot(Ypred,'bs')
    xlabel('q');
    ylabel('S');
    title('Comparing Prediction vs Data');
    legend('Data', 'Model');
end