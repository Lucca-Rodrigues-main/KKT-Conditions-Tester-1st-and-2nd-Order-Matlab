clear all; clc; close all

% Lucca Rodrigues Pinto
% https://github.com/Lucca-Rodrigues-main

% Variaveis
syms x1 x2 x3 real
syms mu1 mu2 mu3 mu4 mu5 mu6 mu7 real
syms y1 y2 y3 real % Corresponde a quantidade de variaveis de x + lambda

% Valores obtidos com KKT de 1 ordem
% [x, lambda, mu]
valvars = [0.2857 1.0179 0.8929 0 0 0 0 0 0 7.0714];

% Funcao objetivo
fo = 2.*x1.^2 + 20.*x2.^2 + 43.*x3.^2 + 12.*x1.*x2 - ...
    16.*x1.*x3 - 56.*x2.*x3 + 8.*x1 + 20.*x2 + 6.*x3;

% Lagrangiano
L = 2.*x1.^2 + 20.*x2.^2 + 43.*x3.^2 + 12.*x1.*x2 - ...
    16.*x1.*x3 - 56.*x2.*x3 + 8.*x1 + 20.*x2 + 6.*x3 ...
    + mu1.*(3.*x1 + 2.*x2 + 5.*x3 - 35) ...
    + mu2.*(-x1 - 2.*x2 - 3.*x3 + 5) ...
    + mu3.*(-x1 + 2.*x2 - 5.*x3 - 3) ...
    + mu4.*(5.*x1 - 3.*x2 + 2.*x3 - 30) ...
    + mu5.*(-x1) ...
    + mu6.*(-x2) ...
    + mu7.*(-x3);

% Condicoes das restricoes de desigualdade
us = [(3.*x1 + 2.*x2 + 5.*x3 - 35); ...
    (-x1 - 2.*x2 - 3.*x3 + 5); ...
    (-x1 + 2.*x2 - 5.*x3 - 3); ...
    (5.*x1 - 3.*x2 + 2.*x3 - 30); ...
    (-x1); (-x2); (-x3)];

% Hessiana do Lagrangiano
variaveis = [x1, x2, x3];
hess_L = hessian(L, variaveis);

hess_L = subs(hess_L, [variaveis,mu1,mu2,mu3,mu4,mu5,mu6,mu7], valvars);

% Condicoes do subespaco M
y = [y1;y2;y3];
for i = 1:length(us)
    M(i) = gradient(us(i), variaveis).' * y == 0;
    M(i) = subs(M(i), variaveis, valvars(1:length(variaveis)));
end
[Sy1,Sy2,Sy3] = vpasolve(M, y.');
Sy = [Sy1,Sy2,Sy3];

expr = y.' * hess_L * y
val = subs(expr, y.', Sy)

if val > 0
    fprintf('Hessiana de L definida positiva em M\n');
elseif val >= 0
    fprintf('Hessiana de L semi-definida positiva em M\n');
elseif val < 0
    fprintf('Hessiana de L definida negativa em M\n');
elseif val <= 0
    fprintf('Hessiana de L semi-definida negativa em M\n');
else
    fprintf('Hessiana de L indefinida em M\n');
end