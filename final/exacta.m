function [exact] = exacta(x,t)

if(t==0)
    exact = ones(size(x));
    exact(x < -1/2) = 0;
    exact(x > 1/2) = 0;
    return;
end

numTerms = 80; 
exact = 1/2*ones(size(x));

for(i=1:numTerms)
    exact=exact+(2/pi)*(sin(i*pi/2)/i)*exp(-i^2*pi^2*t)*cos(i*pi*x);
end