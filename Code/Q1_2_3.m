% NOTE:
% Run Q1_1 before running this file.

Avox = dwis(:,92,65,72);
h=optimset('MaxFunEvals',20000,... 
        'Algorithm','quasi-newton',... 
        'TolX',1e-10,...
        'Display','off',...
        'TolFun',1e-10);

    
    
% Define a starting point for the non­linear fit
startx = [3.5e+00 3e-03 2.5e-01 0 0];
[parameter_hat_final, RESNORM, hessian] = hessianGlobalMin(startx,Avox,bvals,qhat,h);
[diff, predY] = transBallStickSSD(parameter_hat_final,Avox,bvals,qhat);
stdY = sqrt((diff)/(103));
cov = -inv(hessian/(-10000^2));

sigma = sqrt(diag(cov));

twoSigma = [(parameter_hat_final' -sigma), (parameter_hat_final' + sigma)];
laplaceSTD = twoSigma(1:3,:);



%MCMC
[~, paramSamples,~,~] = computeMCMC(parameter_hat_final,Avox,bvals,qhat,diff);
load('../mats/q1_2_1_paramsdist.mat','parameterDists');

for i=1:3
    paramsMean = mean(parameterDists(:,i));
    paramsSTD = std(parameterDists(:,i));
    paramsLen = length(parameterDists(:,i));
    paramsSort = sort(parameterDists(:,i));
    
    % upper/lower limits for 95% confidence
    [lowerLimit, upLimit] = deal(paramsSort(floor(0.025*paramsLen)),paramsSort(ceil(0.975*paramsLen))); 
    
    % Sample STD
    sampleSTD = [paramsMean - paramsSTD, paramsMean + paramsSTD];
    confidence = [lowerLimit, upLimit];
    
    mcmcMean = mean(paramSamples(i,10000:end));
    mcmcSTD = std(paramSamples(i,10000:end));
    mcmcSort = sort(paramSamples(i,10000:end));
    [mcmcLowLim, mcmcUpLim] = deal(mcmcSort(floor(0.025*length(mcmcSort))),mcmcSort(ceil(0.975*length(mcmcSort)))); 
    
    mcmcSampleSTD= [mcmcMean - mcmcSTD, mcmcMean + mcmcSTD];
    mcmcConfidence = [mcmcLowLim, mcmcUpLim];
    
    
    fig = figure;
    plot(mcmcConfidence, [5,5],'g--', 'LineWidth', 0.75);
    ylim([0 6]);
    hold on;
    plot(mcmcSampleSTD, [4,4],'r--', 'LineWidth', 0.75);
    plot(confidence, [3,3],'k:x', 'LineWidth',0.75);
    plot(sampleSTD, [2,2],'b:x', 'LineWidth', 0.75);
    plot(laplaceSTD(i,:), [1,1],'c-o', 'LineWidth', 0.75);
    legend('Confidence limits MCMC','STD MCMC','Confidence limits Parametric','STD Parametric','Laplace','location','northoutside');
    hold off;
    path = sprintf('../Graphs&Maps/Q1_2_3_parameterLimsA%d.png', i);
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

function [acceptanceRate, parameteramples, ratios, newXs] = computeMCMC(x0,Avox,bvals,qhat,startResnorm)
    iters = 100000;
    pertubation = [10000,0.0002,0.01,0.1,0.1];    
    burnin = iters / 10; 
    newParams = ones(5, iters+burnin);
    oldX = x0;
    oldResnorm = startResnorm;

    counter = 0;
    ratios = zeros(1,iters+burnin);
    newXs = zeros(5,iters+burnin);
    for count=1:(iters + burnin)
        newX = normrnd(oldX, pertubation, size(x0));
        newXs(:,count) = newX; 
        [newResnorm, ~] = BallStickSSD(newX,Avox,bvals,qhat);
        acceptanceRatio = ((oldResnorm-newResnorm) / (2*10000^2));
        ratios(count) = acceptanceRatio;
        if (acceptanceRatio > log(rand))
            newParams(:,count) = newX;
            counter = counter + 1;
            oldX = newX;
            oldResnorm = newResnorm;
        else
            newParams(:,count) = oldX;
        end
    end

    acceptanceRate = counter / (iters + burnin);
    parameteramples = newParams;

end


function [finalParams, resnorm, hessian] = hessianGlobalMin(x0, Avox, bvals, qhat,h)

    [x1,x2,x3,x4,x5] = inverse(x0);
    invX = [x1, x2, x3, x4, x5];
    
    % Now run the fitting
    [parameter_hat,RESNORM,~,~,~,hessian] = fminunc('transBallStickSSD', invX, h, Avox, bvals,qhat);
    lowest = RESNORM;

    epochs = 10;
    counter = 0;

    newStart = x0 + 0.0001;
    finalResnorm = ones(1,epochs);
    finalParam = parameter_hat;
    for i=1:epochs
        random = rand(size(newStart));
        [x1,x2,x3,x4,x5] = inverse((newStart) .* random);
        random_start = [x1,x2,x3,x4,x5];
        [parameter_hat,newResnorm,~,~,~,newHessian] = fminunc('transBallStickSSD', random_start, h, Avox, bvals,qhat);
        finalResnorm(i) = newResnorm;
        if (abs(newResnorm - lowest) <= 0.1)
            counter = counter + 1;
        elseif ((lowest - newResnorm) > 0.1)
            lowest = newResnorm;
            counter = 0;
            finalParam = parameter_hat;
            hessian = newHessian;
        end
    end

    finalParams = finalParam;
    resnorm = lowest;
end