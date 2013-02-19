function [exact] = exactb(x,t)

exact = exp(-pi^2*t)*cos(pi*x);