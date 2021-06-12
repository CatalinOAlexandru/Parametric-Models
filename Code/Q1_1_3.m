% NOTE:
% Run Q1_1 before running this file.

Avox = dwis(:,92,65,72);

% Define a starting point for the non-linear fit
startx = [3.5e+00 3e-03 2.5e-01 0 0];
[a,b,c,d,e] = invA(startx);
newstartx = [a,b,c,d,e]; 

% Define various options for the non-linear fitting algorithm

h=optimset('MaxFunEvals',20000,... 
    'Algorithm','quasi-newton',... 
    'TolX',1e-10,...
    'Display','iter',...
    'TolFun',1e-10);

% Now run the fitting
[parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('transBallStickSSD',newstartx,h,Avox,bvals,qhat);
[dif, Ypred] = transBallStickSSD(parameter_hat,Avox,bvals,qhat);

disp(parameter_hat);
disp(RESNORM);

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

function [x1,x2,x3,x4,x5] = invA(x)
    x1 = sqrt(x(1));
    x2 = sqrt(x(2));
    x3 = -log((1/x(3)) - 1);
    x4 = x(4);
    x5 = x(5);
end
