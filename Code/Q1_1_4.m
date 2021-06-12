% NOTE:
% Run Q1_1 before running this file.

Avox = dwis(:,92,65,72);
% Avox = dwis(:,50,65,72);
% Avox = dwis(:,100,50,72);

% Define a starting point for the non-linear fit
startx = [3.5e+00 3e-03 2.5e-01 0 0];
[a,b,c,d,e] = invA(startx);
newstartx = [a,b,c,d,e];

% Define various options for the non-linear fitting algorithm

h=optimset('MaxFunEvals',20000,... 
    'Algorithm','quasi-newton',... 
    'TolX',1e-10,...
    'Display','off',...
    'TolFun',1e-10);

% Now run the fitting
[parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('transBallStickSSD',newstartx,h,Avox,bvals,qhat);
[dif, Ypred] = transBallStickSSD(parameter_hat,Avox,bvals,qhat);

lowest = RESNORM;
epochs = 100;
counter = 0;

allresnorms = zeros(101,1);
allresnorms(1) = RESNORM;

for i=1:epochs
    % creating the random numbers. we add 0.0001 because last 2 values in
    % startx are 0.
    randomNumbers = rand(size(startx));
    [a,b,c,d,e] = invA((startx).*randomNumbers);
    randomStart = [a,b,c,d,e];
    
    % computing the fitting
    [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('transBallStickSSD',randomStart,h,Avox,bvals,qhat);
    allresnorms(i+1) = RESNORM;
    
    % checking RESNORM
    if(abs(lowest - RESNORM) <= 0.1)      
        counter = counter+1;
    elseif((lowest - RESNORM) > 0.1)      
        counter = 0;
        lowest = RESNORM;
    end
    
end


disp(lowest);
disp(counter);
% plot(allresnorms);


function [x1,x2,x3,x4,x5] = invA(x)
    x1 = sqrt(x(1));
    x2 = sqrt(x(2));
    x3 = -log((1/x(3)) - 1);
    x4 = x(4);
    x5 = x(5);
end
