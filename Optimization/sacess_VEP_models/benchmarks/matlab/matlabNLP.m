% Example of mixed integer problem
% This example corresponds to MEIGO toolbox example: ex2.m
% http://gingproc.iim.csic.es/meigo.html

global k1=0.09755988;
global k3=0.0391908;
global k2=0.99*k1;
global k4=0.9*k3;

function [F g]=matlabNLP(X)
F=-x(4);
%Equality constraints
g(1)=x(4)-x(3)+x(2)-x(1)+k4*x(4).*x(6);
g(2)=x(1)-1+k1*x(1).*x(5);
g(3)=x(2)-x(1)+k2*x(2).*x(6);
g(4)=x(3)+x(1)-1+k3*x(3).*x(5);
%Inequality constraint
g(5)=x(5).^0.5+x(6).^0.5;
return
