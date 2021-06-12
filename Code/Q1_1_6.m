% NOTE:
% Run Q1_1 before running this file.

Avox = dwis(:,92,65,72);

% Define various options for the non-linear fitting algorithm
h=optimset('MaxFunEvals',20000,... 
    'Algorithm','quasi-newton',... 
    'TolX',1e-10,...
    'Display','off',...
    'TolFun',1e-10);


[~,~,S0, dmat] = diffTensor(Avox,qhat,bvals);

diff = mean(dmat(:));

[eigenVector,eigenValue] = eig(dmat);

normalisedEiganVal = eigenValue./sum(eigenValue);
f1 = abs(normalisedEiganVal(1) - normalisedEiganVal(2));
f2 = abs(normalisedEiganVal(2) - normalisedEiganVal(3));
f3 = abs(normalisedEiganVal(3) - normalisedEiganVal(1));
f = (1/2)  * (f1 + f2 + f3);

elements = diag(eigenValue);
sortedElements = sort(elements);
eig1 = sortedElements(3);
eig2 = sortedElements(2);
[~, positions] = max(elements);
direction = eigenVector(:,positions);

[phi, theta] = cart2sph(direction(1),direction(2),direction(3));
startx = [S0,diff,f,theta,phi,eig1,eig2];
% ^^^ new start


[a,b,c,d,e] = inverse(startx);
newstartx = [a,b,c,d,e];

% Now run the fitting
[parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('transBallStickSSD',newstartx,h,Avox,bvals,qhat);

lowest = RESNORM;
epochs = 100;
counter = 0;

for i=1:epochs

    % creating the random numbers. we add 0.0001 because last 2 values in
    % startx are 0.
    randomNumbers = rand(size(startx));
    [a,b,c,d,e] = inverse((startx + 0.00001).*randomNumbers);
    randomStart = [a,b,c,d,e];

    % computing the fitting
    [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('transBallStickSSD',randomStart,h,Avox,bvals,qhat);

    % checking RESNORM
    if(abs(RESNORM - lowest) <= 0.1)
        counter = counter+1;
    elseif((lowest - RESNORM) > 0.1)
        counter = 0;
        lowest = RESNORM;
    end
end

disp(lowest);
disp(counter);


function [x1,x2,x3,x4,x5] = inverse(x)
    x1 = sqrt(x(1));
    x2 = sqrt(x(2));
    x3 = -log((1/x(3)) - 1);
    x4 = x(4);
    x5 = x(5);
end


function [x1,x2,x3,x4,x5] = transform(x)
    % Extract the parameters
    % As in transVallStickSSD
    x1 = x(1)^2; 
    x2 = x(2)^2; 
    x3 = 1/(1+exp(-1*x(3)));
    x4 = x(4); 
    x5 = x(5);
end

function [sumRes, resJ, S0, dMat] = diffTensor(Avox,qhat,bvals)
    [qX, qY, qZ] = deal(qhat(1,:), qhat(2,:), qhat(3,:));
 
    y = [ones(1,length(Avox));...
         -bvals.*qX.^2;...
         -2*bvals.*qX.*qY;...
         -2*bvals.*qX.*qZ;...
         -bvals.*qY.^2;...
         -2*bvals.*qY.*qZ;...
         -bvals.*qZ.^ 2];
    yTrans = y';
    
    x = yTrans \ log(Avox);
    S0 = exp(x(1));
    dMat = zeros(3,3);
    dMat(1,:) = [x(2),x(3),x(4)];
    dMat(2,:) = [x(3),x(5),x(6)];
    dMat(3,:) = [x(4),x(6),x(7)];
    
    S = exp(yTrans * x);
    sumRes = sum((Avox-S').^2);
    resJ = S;
end