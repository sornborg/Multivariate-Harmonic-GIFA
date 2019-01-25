function [mngifa, mncva, mnraw, mndft, gam, rho, T2, tau2] = MultiGIFA(matfile, x, y, T, nw, freqs, alph)

%
% MultiGIFA: Version 1.1, written by Andrew Sornborger, 11/25/2007
%
% Reference: Sornborger and Yokoo, NeuroImage
%
% Inputs:
%  filename: file containing imaging data
%      note that the filename and variable name (inside the .mat file) must
%      be the same
%  x: number of pixels in x dimension
%  y: number of pixels in y dimension
%  T: number of frames in imaging data
%  nw: time-bandwidth product
%  freqs: frequencies that harmonic estimates are desired at
%  alph: Type I error bound, set to 0.01 as default
%  datatype: usually 'uint8' or 'uint16' for imaging data
%
% Outputs:
%  mngifa: multivariate, multitaper GIFA mean estimate
%  mncva: multivariate, multitaper CVA mean estimate
%  mnraw: multitaper mean estimate, statistical significance given in T2
%  mndft: discrete Fourier transform estimate
%  gam: eigenvalue from GIFA optimization problem
%  rho: eigenvalue from CVA optimization problem
%  T2: Hotelling's T^2 statistic for multitaper harmonic estimates
%  tau2: statistical significance level for given alph

eval(['load ', matfile, ';']);
eval(['X = ', matfile, ';']);

% number of singular values to keep
ncomps = 200;
disp('SVD to reduce dimensionality');
% first svd to reduce dimensionality of raw data
[psi, a, ev] = pca(X, ncomps);

% number of concentrated tapers
M = 2 * nw - 3;
% generate tapers
[e, v] = dpss(T, nw, M);
% calculate H_m(0)'s
H = sum(e, 1);
% index for even H_m(0)'s (only even H_m(0)'s are non-zero)
idx = [1:2:M];
% calculate constant
H2 = sum(H(idx) * H(idx)');

disp('Tapering SVD eigenvectors');
% taper the svd eigenvectors
tapdat = repmat(e, [1 1 ncomps]);
tapdat = permute(tapdat, [2 1 3]);
for i = 1:2*nw-3
  tapdat(i,:,:) = squeeze(tapdat(i,:,:)) .* a(:,1:ncomps);
end;

disp('Subtracting mean');
% subtract multitaper estimate of mean
mtmeanest = H * squeeze(sum(tapdat, 2))/H2;
a(:,1:ncomps) = a(:,1:ncomps) - repmat(mtmeanest, [T 1]);

% retaper the svd mean subtracted eigenvectors
tapdat = repmat(e, [1 1 ncomps]);
tapdat = permute(tapdat, [2 1 3]);
for i = 1:2*nw-3
  tapdat(i,:,:) = squeeze(tapdat(i,:,:)) .* a(:,1:ncomps);
end;
tapdat = permute(tapdat, [2 1 3]);

% threshold for gifa analysis
tau2 = (M - 1) * finv(1.0 - alph, 2 * (M - 1), 2) * (M - 1);

% loop through frequencies
i = 1;
for freq = freqs
    
  disp(['Frequency: ', num2str(freqs(i))]);
    
  % discrete fourier transform of data
  for j = 1:M
    for k = 1:ncomps
      arr2(j,k) = sum(exp(complex(0,1) * 2 * pi * freq * ([1:T] - 1)/T) .* tapdat(:,j,k)');
    end;
  end;
  
  % discrete fourier transform (untapered) harmonic estimate
  for j = 1:3072
    mndft(i,j) = sum(exp(complex(0,1) * 2 * pi * freq * ([1:T] - 1)/T) .* X(:,j)');
  end;
    
  % svd on J_m
  [psi2, a2, ev2] = pca(arr2, M-1);

  % estimate of mu' (mu in dimensionally reduced space)
  muhat = a2(idx,:)' * H(idx)'/H2;
  % estimate of Z (covariance in dimensionally reduced space)
  K = (a2' - muhat * H) * (a2' - muhat * H)';
  % estimate of signal
  Ks = H2 * muhat * muhat';

  % GIFA operator matrix
  G2 = diag(sqrt(ev2(1:M-1))) * (Ks - tau2 * K) * diag(sqrt(ev2(1:M-1)));

  % diagonalize GIFA operator (remember it's rank one, so only one stat. sig. eigenvector)
  [vG2, dG2] = eig(G2);
  idxvG2 = find(real(diag(dG2)) == (max(real(diag(dG2)))));
  % reconstruct GIFA spatial image
  phiG2 = psi(:,1:ncomps) * psi2 * vG2(:,idxvG2);
  eigenval(1,i) = dG2(idxvG2,idxvG2);

  % solve generalized eigenvalue problem for CVA
  [vC2, dC2] = eig(Ks, K);
  idxvC2 = find(real(diag(dC2)) == (max(real(diag(dC2)))));
  % reconstruct CVA spatial image
  phiC2 = psi(:,1:ncomps) * psi2 * diag(sqrt(ev2(1:M-1)))^-1 * vC2(:,idxvC2);
  eigenval(2,i) = dC2(idxvC2,idxvC2);
  % normalize
  phiC2 = phiC2/sqrt(phiC2' * phiC2);

  gam(i) = real(eigenval(1,i));
  rho(i) = real(eigenval(2,i));
  T2(i) = real(H2 * muhat' * (K\muhat));

  mnraw(i,:) = (psi(:,1:ncomps) * psi2(1:ncomps,:) * diag(sqrt(ev2(1:M-1))) * muhat)';
  mngifa(i,:) = ((mnraw(i,:) * phiG2)' * phiG2)';
  mncva(i,:) = ((mnraw(i,:) * phiC2)' * phiC2)';

  i = i + 1;
   
end;
