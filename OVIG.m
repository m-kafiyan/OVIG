function [ms,Hs,vs,Ws] = OVIG(X,y,epsilon,eps)
[n,D] = size(X);
k = length(unique(y));
m0 = zeros(1,D);
beta0=1;
v0=D;
W0=eye(D);

ms = zeros(k,D);
Hs = zeros(k,D,D);
vs = D * ones(k,1);
Ws = zeros(k,D,D);
for ii=1:k
    for jj=1:D
        Ws(ii,jj,jj) = 1;
        Hs(ii,jj,jj) = 1;
    end
end
nk = zeros(k,1);
ps = zeros(k,1);

xhat        = zeros(k,D);
S           = zeros(k,D,D);
J           = zeros(k,D,D);
for t=1:n
    disp(t);
    yt = y(t);
    xt = X(t,:);
    nk(yt) = nk(yt)+1;
    for i=1:k
        sigg = pinv(vs(i)* (reshape(Ws(i,:,:),D,D)) );
        ps(i) = (nk(yt)/t)*mvnpdf(xt,ms(i,:),(sigg+sigg')./2);
    end
    [~,ypt] = max(ps);
    if t>1  
        xhat(yt,:) = ((t-1)/t)*xhat(yt,:) + (1/t)*xt;
        S(yt,:,:) = reshape(S(yt,:,:),D,D) + (t/(t-1))*((xt-xhat(yt,:))'*(xt-xhat(yt,:)));
    else
        xhat(yt,:) = (1/t)*xt;
        S(yt,:,:) = ((xt-xhat(yt,:))'*(xt-xhat(yt,:)));
    end
    
    l_t = (yt==ypt);
    if l_t==0
        J(yt,:,:) = (xhat(yt,:) - m0)'*(xhat(yt,:) - m0);
        [m,H,W,v] = VIG(t, X, epsilon, m0, beta0, v0, W0, S(yt,:,:), J(yt,:,:), xhat(yt,:));
        
%         m=round(m,10);
%         H=round(H,10);
%         W=(W+W')./2;
%         v=round(v,10);
        
        ms(yt,:)=m;
        Hs(yt,:,:)=H;
        
        vs(yt)=v;
%         Ws(yt,:,:)=W;
        
        [~,err]=cholcov(W,0);
        if err ~= 0
            [VV,DD] = eig(W);
            d=diag(DD);
            d(d<=eps) = eps;
            dc=diag(d);
            W_pd = VV*dc*VV';
            Ws(yt,:,:)=W_pd;
        else
            Ws(yt,:,:)=W;
        end
        
    end
end


end