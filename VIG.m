function [m,H,W,v] = VIG(t, X, epsilon, m0, beta0, v0, W0, S, J, xhat)

i = 1;
[~,D] = size(X);
S = reshape(S,D,D);
J = reshape(J,D,D);
n=t;
m = zeros(D,1);
v=D;
W=eye(D);
H=eye(D);
while(i<100)
    prev_lb = LowerBound(W0, v0, W, v, n, D, beta0, H, S, J);
    m =  (beta0 * m0 + n*xhat) / (beta0 + n);
    H = (beta0 +n) * (v*W);
    v = v0 + n + 1;
    W = pinv(pinv(W0)+(beta0+n)*pinv(H) + S + ((beta0*n)/(beta0+n))*J) ;
    if i>1 && (abs(LowerBound(W0, v0, W, v, n, D, beta0, H, S, J) - prev_lb) < epsilon) 
        break;
    end
    i = i+1;
end
end