function [e,t]=hw2(N,lambda)

%code to solve first-order wave equation
%using D+ in time, D0 in space, with a
%periodic grid function
%t = 0-2pi
%v(j,0)=sin(xj)
%N = num spatial steps
%@Author: Joseph Isaacson

h = 2*pi/(N+1); %spatial step size
k = lambda*h;   %temporal step size

%temporal and spatial partitions:
t = 0:k:2*pi;
x = -h:h:2*pi;

%initial condition:
v(:,1) =f(x);

%exact solution:
Vn(:,1)=f(x);

%error at t0 = 0:
e(1)=0;

for(n=1:length(t)-1)
%     plot(v(2:N+2,n),'k')
%         hold on;
%     plot(Vn(2:N+2,n))
%     drawnow;
%     pause(.1)
%     hold off;
    
    %update x_0-x_N
    for(j=2:N+2)
        D0 = (v(j+1,n)-v(j-1,n))/(2*h);
        v(j,n+1)=v(j,n)+k*D0;
    end
    %update boundary points
    v(1,n+1) = v(N+2,n+1);
    v(N+3,n+1) = v(2,n+1);
    
    %calculate error (program up h norm)!!!!
%     e(n+1)=norm(Vn(:,n)-v(:,n),inf);
    temp =0;
    for(i=1:length(x))
        temp = temp + (Vn(i,n)-v(i,n))^2*h;
    end
    e(n+1)=sqrt(temp);
    Vn(:,n+1)=f(x+n*k);

end
% plot(t,e)
end

function z=f(x)
%       for(i=1:length(x))
          x(x<0) = x(x<0) + 2*pi;
          x(x>2*pi) = x(x>2*pi)-2*pi;
          z = sin(x);
%       end
end



