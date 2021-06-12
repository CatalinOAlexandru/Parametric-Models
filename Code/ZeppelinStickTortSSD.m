function [sumRes, resJ] = ZeppelinStickTortSSD(startx, Avox, bvals,qhat)
    % Extract the parameters
    S0 = startx(1);
    diff = startx(2);
    f = startx(3);
    theta = startx(4);
    phi = startx(5);
    eig1 = startx(6);
    eig2 = (1-f)*eig1;
    
    fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
    fibgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
    Si = exp(-bvals*diff.*(fibgrad.^2));
    Se = exp(-bvals.*(eig2 + (eig1 - eig2).*(fibgrad.^2)));
    S = S0*(f*Si + (1-f)*Se);
    
    % Sum of square differences
    sumRes = sum((Avox - S').^2);
    resJ = S;
end