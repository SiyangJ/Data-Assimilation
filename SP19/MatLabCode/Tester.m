%i = 1;

Ex  = xens-xfm(:,i);
Ey  = yens-yfm;
Bxy = Ex*Ey'/(nens-1);
Byy = cov(yens') + Robsmat;
kgain1 = Bxy*pinv(Byy);

H = h_jacobian(xfm(:,i));
pf = xfcov(:,:,i)+inflmu*trace(xfcov(:,:,i))/ndim*eye(ndim);
pfht = pf*H'; 

Pxy = pfht;
Pyy = H*pfht+Robsmat;

kgain2 = pfht*pinv(H*pfht+Robsmat);