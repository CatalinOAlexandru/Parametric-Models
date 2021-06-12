% NOTE:
% Run Q1_3 before running this file.
std = 0.04;
K = length(meas);



load('../mats/q1_3_1B.mat','bestParamHat','bestResnorm');
ballParams = bestParamHat;
ballRes = bestResnorm;

% load('../mats/q1_3_2A.mat','sumRes');
% diffRes = sumRes;

load('../mats/q1_3_2B.mat','bestParamHat','bestResnorm');
zeppParams = bestParamHat;
zeppRes = bestResnorm;

load('../mats/q1_3_2C.mat','bestParamHat','bestResnorm');
zeppTortParams = bestParamHat;
zeppTortRes = bestResnorm;



aicBall = compAIC(ballRes,5,std);
% aicDiff = compAIC(diffRes,5,std);
aicZepp = compAIC(zeppRes,7,std);
aicZeppTort = compAIC(zeppTortRes,6,std);

disp("AIC");
disp(aicBall);
% disp(aicDiff);
disp(aicZepp);
disp(aicZeppTort);

bicBall = compBIC(K,ballRes,5,std);
% bicDiff = compBIC(K,diffRes,6,std);
bicZepp = compBIC(K,zeppRes,7,std);
bicZeppTort = compBIC(K,zeppTortRes,6,std);

disp("BIC");
disp(bicBall);
% disp(bicDiff);
disp(bicZepp);
disp(bicZeppTort);


% N num param 
% K num data points
% sigma is STD
function aic = compAIC(N, SSD, sigma)
aic = 2*N + SSD/(sigma^2);
end

function bic = compBIC(N, K, SSD, sigma)
bic = N*log(K) + SSD/(sigma^2);
end