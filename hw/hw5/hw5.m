% function [e,t]=hw5(ikd,lambda)
function [efin,hfin]=hw5(ikd,lambda)

index = 0;
for(N =[8,16,32,64])
    index = index+1;
clear v;
clear Vn;
clear e;
%code to solve convection-diffusion equation u_t+a*u_x=eta*u_xx
%a = eta = 1 for our purpose
%periodic grid function
%t = 0-1
%v(j,0)=sin(xj)
%N = num spatial steps
%@Author: Joseph Isaacson

T = 1;          %final time
h = 2*pi/(N+1); %spatial step size
% k = lambda*h^2;   %temporal step size


m = ceil(-log(lambda*h^2)/log(10));
k = 10^(-m);



%temporal and spatial partitions:
t = 0:k:1;
x = -h:h:2*pi;

%initial condition:
v(:,1) =f(x);

%exact solution:
Vn(:,1)=f(x);

%error at t0 = 0:
e(1)=0;

for(n=1:length(t)-1)
%     hold on;
%     plot(x(2:length(x)-1),Vn(2:N+2,n))
%     plot(x(2:length(x)-1),v(2:N+2,n),'b')
%     pause(.1)
    
    %update x_0-x_N
    for(j=2:N+2)
        D0  = (v(j+1,n)-v(j-1,n))/(2*h);
        Dpm = (v(j+1,n)-2*v(j,n)+v(j-1,n))/h^2;
        v(j,n+1)=v(j,n)+k*Dpm-k*D0;
    end
    %update boundary points
    v(1,n+1) = v(N+2,n+1);
    v(N+3,n+1) = v(2,n+1);
    
    %calculate L2 error
    temp =0;
    for(i=1:length(x))
        temp = temp + (Vn(i,n)-v(i,n))^2;
    end
    e(n+1)=sqrt(temp);
    norm(Vn(:,n)-v(:,n));
    
    %find exact solution
    Vn(:,n+1)=f(x-n*k)*exp(-n*k);
    
    efin(index)=e(length(e));
    hfin(index)=h;
end
end

plot(hfin,efin,'o-');

end

function z=f(x)
    T = 2*pi;
    z = sin(mod(x+T,T));
end



