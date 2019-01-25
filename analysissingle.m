% the code below loads in the .mat file 'orientrot', a file of 1800 frames of 64 by
%  48 pixel images. using a time-bandwidth product nw = 20, it obtains
%  multivariate, multitaper estimates of the harmonic content at wavenumbers
%  50 and 100. the parameter alpha, used in the GIFA estimate is set to
%  0.01. for further information, type 'help MultiGIFA'.
[mng, mnc, mnr, mnf, gam, rho, T2, tau2] = MultiGIFA('orientrotsingle', 64, 48, 1800, 20, [50, 100], 0.01);

% plot results.
figure;
subplot(4,2,1);
pcolor(reshape(real(mng(1,:)), [64 48])); shading interp; colormap gray;
axis equal; axis tight;
subplot(4,2,2);
pcolor(reshape(imag(mng(1,:)), [64 48])); shading interp; colormap gray;
axis equal; axis tight;
subplot(4,2,3);
pcolor(reshape(real(mnc(1,:)), [64 48])); shading interp; colormap gray;
axis equal; axis tight;
subplot(4,2,4);
pcolor(reshape(imag(mnc(1,:)), [64 48])); shading interp; colormap gray;
axis equal; axis tight;
subplot(4,2,5);
pcolor(reshape(real(mnr(1,:)), [64 48])); shading interp; colormap gray;
axis equal; axis tight;
subplot(4,2,6);
pcolor(reshape(imag(mnr(1,:)), [64 48])); shading interp; colormap gray;
axis equal; axis tight;
subplot(4,2,7);
pcolor(reshape(real(mnf(1,:)), [64 48])); shading interp; colormap gray;
axis equal; axis tight;
subplot(4,2,8);
pcolor(reshape(imag(mnf(1,:)), [64 48])); shading interp; colormap gray;
axis equal; axis tight;

figure;
subplot(4,2,1);
pcolor(reshape(real(mng(2,:)), [64 48])); shading interp; colormap gray;
axis equal; axis tight;
subplot(4,2,2);
pcolor(reshape(imag(mng(2,:)), [64 48])); shading interp; colormap gray;
axis equal; axis tight;
subplot(4,2,3);
pcolor(reshape(real(mnc(2,:)), [64 48])); shading interp; colormap gray;
axis equal; axis tight;
subplot(4,2,4);
pcolor(reshape(imag(mnc(2,:)), [64 48])); shading interp; colormap gray;
axis equal; axis tight;
subplot(4,2,5);
pcolor(reshape(real(mnr(2,:)), [64 48])); shading interp; colormap gray;
axis equal; axis tight;
subplot(4,2,6);
pcolor(reshape(imag(mnr(2,:)), [64 48])); shading interp; colormap gray;
axis equal; axis tight;
subplot(4,2,7);
pcolor(reshape(real(mnf(2,:)), [64 48])); shading interp; colormap gray;
axis equal; axis tight;
subplot(4,2,8);
pcolor(reshape(imag(mnf(2,:)), [64 48])); shading interp; colormap gray;
axis equal; axis tight;
