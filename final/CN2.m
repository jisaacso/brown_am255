function err = CN2(h,mu)
%h = spatial step size
%mu = k/h^2;

x = -1:h:1;
N = length(x);

% k = h^2*mu;
k = 1 / (2 / h);
t = [0:k:1];
T = length(t);
alpha = k/(2*h^2);

vn = zeros(N,1);
vn = exactb(x,0);
vn = vn';
v0=vn(1:(end-1));

Q2 = diag(alpha*ones(N-2,1),1)+diag(alpha*ones(N-2,1),-1)+diag(ones(N-1,1)-2*alpha*ones(N-1,1));
Q2(1,N-1)=alpha;
Q2(N-1,1)=alpha;
Q1 = -Q2;
Q1 = Q1+diag(2*ones(N-1,1));


Q = Q1\Q2;


for(i = 1:(T-1))
   vn = Q*v0;
   v0 = vn;
end

err = norm((exactb(x(1:(end-1)),t(T))'-vn)*sqrt(h));