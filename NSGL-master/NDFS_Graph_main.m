clc
clear all;
warning off;
load(' ');data = ' ';v = size(X,2);c = ;

XX = DataConcatenate(X);
[n,d] = size(XX);
XX = XX';
% parameters 
param.beta = 1e+4;
param.gamma = 1e+4;
param.lambda = 1e-4;
param.v = v;
param.t = 2;
param.k = 8;
param.n = n;
param.d = d;
param.c = c;
param.NITER = 100;
l = 100;
a = 1;
p = 1;

[F,W,Wi,S,L,M,Wv,iter] = NDFS_Graph(XX,X,param,alpha);
d = size(Wi,1);
[Wi_des, Wi_index] = sort(Wi,'descend');
Wi_idx = Wi_index(1:l,:);
Wi_l = Wi_des(1:l,:);
Wi_identify = zeros(l,d);
for i = 1:l
    Wi_identify(i,Wi_idx(i)) = 1;
end
Xw = Wi_identify*XX;
for m = 1:10
    [y] = litekmeans(Xw', param.k);
    result = ClusteringMeasure(Y,y);
    fprintf('dataset = %s, k = %d, beta = %d,gamma = %d, lambda = %d, ',...
        data,param.k,alpha,param.beta,param.gamma,param.lambda);
    fprintf('\n');
    disp(['Best. ACC: ',num2str(result(1))]);
    disp(['Best. NMI: ',num2str(result(2))]);
    disp(['Best. Purity: ',num2str(result(3))]);
    fprintf('\n');
    Fin_(a,(1:3)) = result;
    a = a+1;
    Fin_result(m,(1:3)) = result;
end
result1 = sum(Fin_result);
result2 = result1/10;