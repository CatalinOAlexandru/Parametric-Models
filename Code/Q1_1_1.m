% NOTE:
% Run Q1_1 before running this file.

% Question 1.1.1 as seen in the Lecture 2a

% Extract the set of 108 measurements from one single voxel
Avox = dwis(:,92,65,72);

G = [ones(1,108);
     -bvals.*qhat(1,:).^2;
     -2*bvals.*qhat(1,:).*qhat(2,:);
     -2*bvals.*qhat(1,:).*qhat(3,:);
     -bvals.*qhat(2,:).^2;
     -2*bvals.*qhat(2,:).*qhat(3,:);
     -bvals.*qhat(3,:).^2]';

x = G\log(Avox);

% difusion tensor
D = [[x(2) x(3) x(4)]; [x(3) x(5) x(6)]; [x(4) x(6) x(7)]];

% mean difusibity
md = trace(D)/3;

% eigen vectors (a) and eigen values (b)
[a,b] = eig(D);

dt_map72 = zeros(7, 145, 174);
FA = zeros(145,174);

Ginv = pinv(G);

for i=1:145
    for j=1:174
        A = dwis(:,i,j,72);
        if(min(A)>0)
%           dt_map72(:,i,j) = G\log(A);
            dt_map72(:,i,j) = Ginv * log(A);  
            x = G\log(Avox);
            D = [[x(2) x(3) x(4)]; [x(3) x(5) x(6)]; [x(4) x(6) x(7)]];
            FA(i,j) = sqrt(1/2 * (3 - 1/trace((D/trace(D))^2)));
        end
    end
end


save('../mats/q1_1_1.mat','FA');