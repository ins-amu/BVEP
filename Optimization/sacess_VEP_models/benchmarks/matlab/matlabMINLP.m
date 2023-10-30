% Example of mixed integer problem
% This example corresponds to MEIGO toolbox example: ex3.m
% http://gingproc.iim.csic.es/meigo.html

function [F g]=matlabMINLP(x)
F = x(2)^2 + x(3)^2 + 2.0*x(1)^2 + x(4)^2 - 5.0*x(2) - 5.0*x(3) - 21.0*x(1) + 7.0*x(4);
g(1) = x(2)^2 + x(3)^2 + x(1)^2 + x(4)^2 + x(2) - x(3) + x(1) - x(4);
g(2) = x(2)^2 + 2.0*x(3)^2 + x(1)^2 + 2.0*x(4)^2 - x(2) - x(4);
g(3) = 2.0*x(2)^2 + x(3)^2 + x(1)^2 + 2.0*x(2) - x(3) - x(4);
return
