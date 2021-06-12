format long;

load('../mats/q1_3_1B.mat','bestParamHat','bestResnorm');

Avox = meas;
qhat = grad_dirs;
sigma = 0.04;

length_params = length(bestParamHat);
length_input = length(Avox);
fisher = zeros(length_params,length_params);
derivs = getDerivatives(bestParamHat,bvals,qhat);
for i=1:length_input
    fisher = fisher + derivs(i,:)'*derivs(i,:);
end
fisher = fisher./(sigma^2);
fisher = fisher./(bestParamHat'*bestParamHat);

disp(fisher);
save('../mats/q1_4_1.mat','fisher');

function [derivs] = getDerivatives(startx, bvals,qhat)
    % Extract the parameters
    S0 = startx(1);
    diff = startx(2);
    f = startx(3);
    theta = startx(4);
    phi = startx(5);
    
    fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
    fibgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');    
    Si = exp(-bvals * diff .* (fibgrad.^2));
    Se = exp(-bvals*diff);
    deriv_F = S0*(Si-Se);
    
    intermediate = S0*f*Si .* (-(2*bvals*diff) .* fibgrad);
    ddTheta = repmat([cos(phi)*cos(theta),sin(phi)*cos(theta),-sin(theta)], [length(qhat) 1]);
    ddPhi = repmat([-sin(phi)*sin(theta),cos(phi)*sin(theta),0], [length(qhat) 1]);
    derPhi = intermediate.*sum(ddPhi'.*qhat);
    derTheta = intermediate.*sum(ddTheta'.*qhat);
    
    derDiff = S0*(f*Si.*(-bvals*diff.*(fibgrad.^2))) + ((1-f)*Se.*(-bvals));
    derS0 = f*Si + (1-f)*Se;
    derivs = [derS0;derDiff;deriv_F;derTheta;derPhi]';
end