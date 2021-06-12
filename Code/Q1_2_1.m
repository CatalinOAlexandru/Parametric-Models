% NOTE:
% Run Q1_1 before running this file.

Avox = dwis(:,92,65,72);
% Avox = dwis(:,50,65,72);
h=optimset('MaxFunEvals',20000,... 
        'Algorithm','quasi-newton',... 
        'TolX',1e-10,...
        'Display','off',...
        'TolFun',1e-10);

% Define a starting point for the non­linear fit
% New better start
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



% From Q1_1_4
[a,b,c,d,e] = inverse(startx);
newstartx = [a,b,c,d,e];

% Now run the fitting
[parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('transBallStickSSD',newstartx,h,Avox,bvals,qhat);

parameter_hat_final = parameter_hat;
lowestResnorm = RESNORM; % global min

epochs = 5;
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
    if(abs(RESNORM - lowestResnorm) <= 0.1)
        counter = counter+1;
    elseif((lowestResnorm - RESNORM) > 0.1)
        counter = 0;
        lowestResnorm = RESNORM;
        parameter_hat_final = parameter_hat;
    end
end


[diff, predY] = transBallStickSSD(parameter_hat_final,Avox,bvals,qhat);
stdY = sqrt(diff/103);

iters = 1000;
parameterDists = zeros(iters,5);

for i=1:iters
    normPred = predY + normrnd(0,stdY, size(predY));
    [newParamHat,res,exit,out] = fminunc('transBallStickSSD', parameter_hat_final, h, normPred', bvals,qhat);
    [p1,p2,p3,p4,p5] = transform(newParamHat);
    parameterDists(i,:) = deal([p1,p2,p3,p4,p5]);
end


save('../mats/q1_2_1_paramsdist.mat','parameterDists');


for i=1:size(parameterDists(1,:),2)
    paramsMean = mean(parameterDists(:,i));
    paramsSTD = std(parameterDists(:,i));
    paramsLen = length(parameterDists(:,i));
    paramsSort = sort(parameterDists(:,i));
    
    % upper/lower limits for 95% confidence
    [lowerLimit, upLimit] = deal(paramsSort(floor(0.025*paramsLen)),paramsSort(ceil(0.975*paramsLen))); 
    
    % Sample STD
    sampleSTD = [paramsMean - paramsSTD, paramsMean + paramsSTD];
    confidence = [lowerLimit, upLimit];
    
    fig = figure;
    hist1 = histfit(paramsSort,50);
    boundY = ylim;
    hold on;
    set(hist1(1),'facecolor','none'); set(hist1(2),'color','m');
    plot(confidence, 0.5*[boundY(2) ,boundY(2)],'b-', 'LineWidth',0.75);
    plot(sampleSTD, 0.6*[boundY(2) ,boundY(2)],'r-', 'LineWidth', 0.75);
    legend('Likelihood distribution','Likelihood Fit','95% Range','2 Sigma','location','northoutside');
    
    path = sprintf('../Graphs&Maps/Q1_2_1_parameterA%d.png', i);
    saveas(fig,path);
end




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
