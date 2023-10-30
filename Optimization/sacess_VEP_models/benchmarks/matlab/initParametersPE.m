% In this function, you can define like a global variable Matlab data
% If you define global variables, remember declaring these variables in the corresponding 
% matlab objective function.

function [e]=initParametersPE(E)
global texp  yexp

texp= [0 1230 3060 4920 7800 10680 15030 22620 36420];

% Distribution of species concentration
% y(1) y(2) y(3) y(4) y(5)
yexp= [ 100.0 0.0 0.0 0.0 0.0
        88.35 7.3 2.3 0.4 1.75
        76.4 15.6 4.5 0.7 2.8
        65.1 23.1 5.3 1.1 5.8
        50.4 32.9 6.0 1.5 9.3
        37.5 42.7 6.0 1.9 12.0
        25.9 49.1 5.9 2.2 17.0
        14.0 57.4 5.1 2.6 21.0
        4.5 63.1 3.8 2.9 25.7 ];

e=0;
return
