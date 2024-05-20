clc;
clear all;
close all;
load vpdataOLS;

%% Finding the initial estimates A_in and B_in using linear regression

N = length(psat);
y = log(psat); % Divided by 100 as converting to bar
temp = temp + 273;
x = 1./(temp+273);

covxy = cov(x,y);
meanx = mean(x);
meany = mean(y);
B_in = covxy(1,2)/covxy(1,1); % This is initial guess of B
A_in = meany - B_in*meanx; % This is initial guess of A

x0(1) = A_in;
x0(2) = -B_in;
x0(3) = 273;

options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','Display','iter');
x_1 = lsqnonlin(@(x)fun(x,psat,temp),x0,[],[],options); % This is the non linear optimization function fun

x_2 = fmincon(@(x)fun1(x,psat,temp),x0, [], []); % This is the function fmincon which optimizes scalar returned by fun1

%% Printing and plotting
clc;

disp("The values obtained from fmincon are as follows: "); disp("");
fprintf("A: %0.4f\n", x_1(1));
fprintf("B: %0.4f\n", x_1(2));
fprintf("C: %0.4f\n", x_1(3));

disp("The values obtained from lsqnonlin are as follows: "); disp("");
fprintf("A: %0.4f\n", x_2(1));
fprintf("B: %0.4f\n", x_2(2));
fprintf("C: %0.4f\n", x_2(3));

scatter(temp, x_1(1) - x_1(2)./(temp + x_1(3)), "*", "red");
hold on;
scatter(temp, x_2(1) - x_2(2)./(temp + x_2(3)), "+", "blue");
legend("lsqnonlin", "fmincon");
xlabel("Temperature in Kelvin");
ylabel("ln(p^{sat}) in bar");
title("ln(p^{sat}) vs Temperature as per fmincon, lsqnonlin")
hold off;

function [res] = fun(x,y,T)
N = length(T);
res = zeros(N,1);
for i = 1:N
    res(i) = y(i) - exp(x(1) - x(2)/(T(i) + x(3)));
end
end

function res1 = fun1(x,y,T)
N = length(T);
res1 = 0;
for i = 1:N
    res1 = res1 + (y(i) - exp(x(1) - x(2)/(T(i) + x(3))))^2;
end
end