% clear;clc;
% a=readtable('abalone.csv');
% [~, ~, y] = unique(a{:,1});
% X=a{:,2:end};

clear;clc;
[X,y] = wine_dataset;
X=X';
y=y';
idx=randperm(length(y));
X=X(idx,:);
y=y(idx,:);
[~,y] = max(y,[],2);

X = (X - repmat(mean(X),size(X,1),1) ) ./ repmat(std(X),size(X,1),1) ;
Xtest = X(1:20,:);
ytest = y(1:20);
Xtrain = X(20:end,:);
ytrain = y(20:end);

eps = 1e-20;
epsilon = 1e-20;
[ms,Hs,vs,Ws] = OVIG(Xtrain,ytrain,epsilon,eps);

k=max(y);
[n,D] = size(Xtest);
ypt=zeros(n,1);
for ind=1:n
    for i=1:k
        ps(i) = mvnpdf(Xtest(ind,:),ms(i,:),pinv(vs(i)* (reshape(Ws(i,:,:),D,D)) ));
    end

    [~,ypt(ind)] = max(ps);
end

fprintf('error = ')
disp(sum(ytest~=ypt) / n);