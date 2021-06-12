% NOTE:
% Run Q1_1 before running this file.

Avox = dwis(:,92,65,72);

% Define various options for the non-linear fitting algorithm
h=optimset('MaxFunEvals',20000,... 
    'Algorithm','quasi-newton',... 
    'TolX',1e-10,...
    'Display','off',...
    'TolFun',1e-10);

% Define a starting point for the non-linear fit
startx = [3.5e+00 3e-03 2.5e-01 0 0];

[images, voxelX, voxelY, slices] = size(dwis);
allParameters = zeros(slices,voxelX,5);
allResnorms = zeros(slices,voxelX,5);

tic;

for j=1:voxelX
   for k=1:voxelY
       
       fprintf('J: %d/%d | K: %d/%d \n',j,voxelX,k,voxelY);
       
       voxel = dwis(:,j,k,72);
       
       [a,b,c,d,e] = invA(startx);
       newstartx = [a,b,c,d,e];
       
       try
           % Now run the fitting
           [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('transBallStickSSD',newstartx,h,voxel,bvals,qhat);
       catch
           disp('found error');
           continue
       end
       
       parameter_hat_final = parameter_hat;
       lowest = RESNORM; % global min
       
       epochs = 5;
       counter = 0;

       for i=1:epochs
            
            % creating the random numbers. we add 0.0001 because last 2 values in
            % startx are 0.
            randomNumbers = rand(size(startx));
            [a,b,c,d,e] = invA((startx + 0.00001).*randomNumbers);
            randomStart = [a,b,c,d,e];

            % computing the fitting
            [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('transBallStickSSD',randomStart,h,voxel,bvals,qhat);

            % checking RESNORM
            if(abs(RESNORM - lowest) <= 0.1)
                counter = counter+1;
            elseif((lowest - RESNORM) > 0.1)
                counter = 0;
                lowest = RESNORM;
                parameter_hat_final = parameter_hat;
            end
       end
       
       [x1,x2,x3,x4,x5] = transform(parameter_hat_final);
       allParameters(j,i,:) = [x1,x2,x3,x4,x5];
       allResnorms(j,i,:) = lowest;

   end
end

toc;


% Saving mats
% save('../mats/q1_1_5_params.mat','allParameters');
% save('../mats/q1_1_5_resnorms.mat','allResnorms');



function [x1,x2,x3,x4,x5] = transform(x)
    % Extract the parameters
    % As in transVallStickSSD
    x1 = x(1)^2; 
    x2 = x(2)^2; 
    x3 = 1/(1+exp(-1*x(3)));
    x4 = x(4); 
    x5 = x(5);
end


function [x1,x2,x3,x4,x5] = invA(x)
    x1 = sqrt(x(1));
    x2 = sqrt(x(2));
    x3 = -log((1/x(3)) - 1);
    x4 = x(4);
    x5 = x(5);
end
