% Example of parameter estimation problem
% This example corresponds to MEIGO toolbox example: ex5.m
% http://gingproc.iim.csic.es/meigo.html



function [F R]=matlabPE(X)
	global texp yexp

	[tout,yout] = ode15s(@ex5_dynamics,texp,[100 0 0 0 0],[],X);
	R=(yout-yexp);
	R=reshape(R,numel(R),1);
	F = sum(sum((yout-yexp).^2));
        save X
        save R
	% Vectorial outputs always need to be loaded in rows.
	% You need to return F and R
        R=R';
return
%***************************************************
%Function of the dynamic system
function dy=ex5_dynamics(t,y,p)
	dy=zeros(5,1); %Initialize the state variables
	dy(1)=-(p(1)+p(2))*y(1);
	dy(2)=p(1)*y(1);
	dy(3)=p(2)*y(1)-(p(3)+p(4))*y(3)+p(5)*y(5);
	dy(4)=p(3)*y(3);
	dy(5)=p(4)*y(3)-p(5)*y(5);
return
%***************************************************
