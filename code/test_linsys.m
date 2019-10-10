clearvars; close all; clc;

% generate data
q = 8; dt = 2*pi/q; r = 2; s = 1e-2; n = 100; rho = 10;

x0 = [1 0.1]';
A = [1 -2; 1 -1];
eA = expm(dt*A);

CZ = zeros(r,q); CZ(:,1) = x0;
for j = 2:q
    CZ(:,j) = eA*CZ(:,j-1);
end

% monte carlo simulation
el = zeros(r,n);
for i = 1:n
    NU = sqrt(s) .* randn( size(CZ) );
    Z = CZ + NU;

    X = Z(:,1:end-1); 
    Y = Z(:,2:end);

    cdmd = CDMD_ADMM(X,Y,'r',r,'rho',rho);
    [EV,EL,A,B,U] = cdmd.compCDMD();

     el(:,i) = sort( log(diag(EL)) / dt );
end
mel = mean(el,2);

% eigenvalues plut
figure;
plot( 0, 1, '^k', 'markers', 30 );
hold on; plot( real(mel(end)), imag(mel(end)), '.k', 'markers', 30 );
plot( real(el(end,:)), imag(el(end,:)), '.', 'color', [0.5,0.5,0.5], 'markers', 20 );
legend({'ground-truth','numerical mean'});
th = 0.05; axis([-th th 1-th 1+th]);