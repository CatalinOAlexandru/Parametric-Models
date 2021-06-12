% NOTE:
% Run Q1_3 before running this file.

Avox = meas;
qhat = grad_dirs;

% Define various options for the non-linear fitting algorithm
h=optimset('MaxFunEvals',20000,... 
    'Algorithm','quasi-newton',... 
    'TolX',1e-10,...
    'Display','off',...
    'TolFun',1e-10);


% Define a starting point for the non­linear fit
% New better start
[sumRes, resJ, S0, dMat] = diffTensor(Avox,qhat,bvals);

disp(sumRes);
plotting(resJ, Avox)


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
