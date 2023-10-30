% Example of constrained problem.
% This example corresponds to MEIGO toolbox example: ex2.m
% http://gingproc.iim.csic.es/meigo.html

function [y g R] = matlabCost(x)
  g = [];
  R = []; 
  % specific function code begins
  y=-x(1)-x(2)
  g(1)=x(2)-2*x(1).^4+8*x(1).^3-8*x(1).^2;
  g(2)=x(2)-4*x(1).^4+32*x(1).^3-88*x(1).^2+96*x(1); 
  % specific function code finalizes
  % always return: 
  %   (1)  y = fitness of the function.
  %   (2)  R = residual of the function. 
  %   (3)  g = penalty of the function.
end


