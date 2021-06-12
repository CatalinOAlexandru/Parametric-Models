% NOTE:
% Run Q1_1_5 before running this file.

load('../mats/q1_1_5_params.mat');
load('../mats/q1_1_5_resnorms.mat');

S0 = allParameters(:,:,1);
d = allParameters(:,:,2);
f = allParameters(:,:,3);

fig = figure;
imshow(S0,[1, inf]);
saveas(fig,'../Graphs&Maps/Q1_1_5_S0.png');
clf(fig);

imshow(d,[1e-04,5e-03]);
saveas(fig,'../Graphs&Maps/Q1_1_5_d.png');
clf(fig);

imshow(f, [min(f(:)),max(f(:))]);
saveas(fig,'../Graphs&Maps/Q1_1_5_f.png');
clf(fig);

imshow(allResnorms, [1e04,max(allResnorms(:))]);
saveas(fig,'../Graphs&Maps/Q1_1_5_resnorms.png');
clf(fig);

theta = allParameters(:,:,4);
phi = allParameters(:,:,5);
xDir = f.*cos(phi).*sin(theta);
yDir = f.*sin(phi).*sin(theta);    

fig = quiver(xDir,yDir);
saveas(fig,'../Graphs&Maps/Q_1_5_quiver.png');