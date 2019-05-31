function [S]  = UpdateDiffusionMatrixSquareRoot(A,X,NumChannels)

N = length(X);
O_N = ones(N,1);
I_N = eye(N,N);

%construct D from A and X, according to equation (11) in Schmerl and McDonnell, Physical Review E
D = ((A*X)*O_N') .* I_N - A.*(O_N*X') - A' .* (X*O_N');
D = D/NumChannels;

%find the matrix square roots; this method assumes D is symmetric, which it is here
[u,s,v] = svd(D);
S = u*sqrt(s)*v';



