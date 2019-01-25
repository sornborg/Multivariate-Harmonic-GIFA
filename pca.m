function [psi, a, evals] = pca(f, N, type)

% [psi, a, evals]=pca(f, N, type)
%
% Version 1.1, for inclusion with SOARS version 1.1 package. In the
% previous version of this code, one of the sets of eigenvectors (psi) was
% not correctly normalized.
%
% performs KL analysis on snaps 
% Inputs:
%   f: data a T x M matrix (smallest dimension first)
%   N: the number of eigenfunctions to calculate
% type is 'float32' [default], 'int16' or other
%
% returns:
%   psi: matrix whose columns are spatial eigenvectors
%   a:   matrix whose columns are temporal eigenvectors
%   evals:  vector of eigenvalues in descending order

% assume floats unless told otherwise
if nargin < 5
  type = 'float32';
end

% find number of pixels (assumed that data is in time x pixel format)
P = size(f,2);

% calculate time-time covariance matrix
cor = f * f';
% find eigenvectors and eigenvalues of time-time covariance matrix
[a, lam] = eig(cor);
% sort into ascending order
[evals, ind] = sort(diag(lam));
% change to descending order
evals = sqrt(flipud(evals));
% save only the eigenvalues (they are on the diagonal)
lam = diag(evals);
% find singular values (square root of eigenvalues)
isig = diag(evals.^(-.5));
% rearrange time eigenvectors  to descending order
a = a(:,flipud(ind));
% only keep N of them
a = a(:,1:N);
% find first N left eigenvectors (images) by inverting decomposition
psi = f' * a * isig(1:N,1:N)^2;
