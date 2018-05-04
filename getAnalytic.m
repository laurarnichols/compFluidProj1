function [ Tact ] = getAnalytic( X, Y, T )
% Supposed to give analytic solution but doesn't work

Tact = zeros(size(T));

for n = 1:1000
    add1 = sin(n*pi*X/0.5)*(500*cosh(n*           pi*Y/0.5) - ...
        500*coth(n*pi)*sinh(n*pi*Y/0.5));
    %input('1');
    add2 = 200/(n*pi*sinh(n*pi))*(1-cos(n*pi))*sinh(n*pi*X/0.5)*sin(n*pi*Y/0.5);
    %input('2');
    add3 = 200*0.5/(n*pi)^2*(1-cos(n*pi))* ...
        (sinh(n*pi*X/0.5) - tanh(n*pi)*cosh(n*pi*X/0.5))* ...
        sin(n*pi*Y/0.5);
    %input('3');
    add4 = 5*0.5/((n*pi)^2*cosh(n*pi))*(cos(n*pi) - 1)* ...
           sinh(n*pi*Y/0.5)*sin(n*pi*X/0.5);
    %input('4');
    Tact = Tact + add1 + add2 + add3 + add4; 
end

