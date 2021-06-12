

load('../data/data'); 
dwis=double(dwis); 
dwis=permute(dwis,[4,1,2,3]);
qhat = load('../data/bvecs');
bvals = 1000*sum(qhat.*qhat);

% imshow(flipud(squeeze(dwis(1,:,:,72))'), []);