function [F,W,Wi,S,L,M,Wv,iter] = NDFS_Graph(XX,X,param,alpha)
%% ===================== Parameters =====================
beta = param.beta;
gamma = param.gamma;
lambda = param.lambda;
NITER = param.NITER;
v = param.v;
k = param.k;
n = param.n;
d = param.d;
t = param.t;
c = param.c;
%% =====================   Normalization =====================
% for i = 1 :v
%     for  j = 1:n
%         X{i}(:,j) = ( X{i}(:,j) - mean( X{i}(:,j) ) ) / std( X{i}(:,j) ) ;
%         Xx{i} = X{i}*X{i}';
%     end
% end
%% ===================== initialize =====================
XtX = XX*XX';
Wv= 1/v*ones(1,v);
W = zeros(d,c);
Wi = sqrt(sum(W.*W,2)+eps);
dd = 0.5./Wi;
W_D = diag(dd);
SUM = zeros(n);
for i = 1:v
    distX_initial(:,:,i) =  L2_distance_1( X{i}',X{i}' ) ;                  %initialize X
    SUM = SUM+distX_initial(:,:,i);
end
distX = 1/v*SUM;
[distXs, idx] = sort(distX,2);

%initialize S
S = zeros(n);
rr = zeros(n,1);
for i = 1:n
    di = distXs(i,2:k+2);
    rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
    id = idx(i,2:k+2);
    S(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end;
S = (S+S')/2;
D = diag(sum(S));
L = D-S;
[F,~,~] = eig1(L,c);
% [F0, ~, evs]=eig1(L0, n, 0);
% F = F0(:,2:(c+1));

%% ===================== updating =====================
for iter = 1:NITER
    % update weighted_distX
    SUM = zeros(n,n);
    for i = 1 : v
        if iter == 1
            distX_updated = distX_initial;
        end
        % update X
        distX_updated(:,:,i) = Wv(i)*distX_updated(:,:,i) ;
        SUM = SUM + distX_updated(:,:,i);
    end
    distX = SUM;
    
    %update S
    distf = L2_distance_1(F',F');
    S = zeros(n);
    for i=1:n
        idxa0 = idx(i,2:k+1);
        dfi = distf(i,idxa0);
        dxi = distX(i,idxa0);
        ad = -(dxi+lambda*dfi)/(2*alpha);
        S(i,idxa0) = EProjSimplex_new(ad);
    end;
    S = (S+S')/2;
    D = diag(sum(S));
    L = D-S;
    
    sum_M = zeros(n);
    for i = 1:v                                                         %update W
        G = inv(XtX+W_D);
        W = G*XX*F;
        Wi = sqrt(sum(W.*W,2)+eps);
        dd = 0.5./Wi;
        W_D = diag(dd);
        %         clear Wi
        M = lambda*L+beta*(eye(n)-XX'*G*XX);
        clear G
        sum_M = sum_M+Wv(i)*M;
    end
    M = (sum_M+sum_M')/2;
    
    F = F.*(gamma*F + eps)./(M*F + gamma*F*F'*F + eps);                 %update F
    F = F*diag(sqrt(1./(diag(F'*F)+eps)));
    
    for i = 1:v
        h(i) = trace(X{i}'*S*X{i});
        H(i) = (1/h(i)).^(1/(t-1));
    end
    H1 = sum(H);
    
    for i = 1:v
        Wv(i) = H(i)./H1;
    end
    Obj = 0;
    for i = 1:v
        obj = Wv(i)*trace(X{i}'*S*X{i});
        Obj = Obj+obj;
    end
    OBJ(iter) = Obj+alpha*norm(S,'fro').^2+trace(F'*M*F)+gamma/4*norm(F'*F-eye(size(F,2)),'fro')^2;
    if iter == 1
        err = 0;
    else
        err = OBJ(iter-1)-OBJ(iter);
    end
    
    fprintf('iteration =  %d:  obj: %.4f; err: %.4f  \n', ...
        iter, OBJ(iter), err);
    if (abs(err))<1e-2
        if iter > 2
            break;
        end
    end
end




















