function err = euler(h,mu,idx)
%h = spatial step size
%mu = k/h^2;

x = -1:h:1;
N = length(x);

k = h^2*mu;
% k = (.0025)*(1/4)^(idx-1);
t = [0:k:1,1];
% t = 0:k:1;
% mu = k/(h^2);
T = length(t);

v = zeros(N,T);
v(:,1) = exactb(x,0);

Q = zeros(N,N);
Q = Q+diag(mu*ones(N-1,1),1)+diag(mu*ones(N-1,1),-1)+diag(ones(N,1)-2*mu*ones(N,1));
Q(1,N-1)=mu;
Q(N,2)=mu;

for(i = 1:T-1)
      v(:,i+1) = Q*v(:,i);
end

err = norm((exactb(x,t(T))'-v(:,T))*sqrt(h));


