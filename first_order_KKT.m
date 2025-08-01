clear all; clc; close all

% Lucca Rodrigues Pinto
% https://github.com/Lucca-Rodrigues-main

% Variaveis
syms x1 x2 x3 real
syms mu1 mu2 mu3 mu4 mu5 mu6 mu7 real

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

% Restricoes de igualdade
lambdas = [];

% Gradiente do Lagrangiano
variaveis = [x1, x2, x3];
va = 'x1, x2, x3';
grad_L = gradient(L, variaveis);

% Gerar todas as combinacoes de equacoes ativas e inativas
mu_comb = dec2bin(0:(2^length(us))-1) - '0';
k = 1;

% Iterar sobre cada combinação
for i = 1:size(mu_comb, 1)
    disp(i);
    % Indices das equacoes ativas
    ativas = find(mu_comb(i, :) == 1);
    % Indices das equacoes inativas
    inativas = find(mu_comb(i, :) == 0);
    
    % Usa solve para resolver cada caso
    if length(inativas) > 0 && length(ativas) > 0
        temp = sprintf('mu%d,',inativas);
        temp2 = sprintf('%d,',zeros(1,length(inativas)));
        temp3 = sprintf('mu%d,',ativas);
        
        try
            eval(['sol = struct2cell(vpasolve(subs([grad_L ; us(ativas)], [' temp(1:end-1) '], [' temp2(1:end-1) '])==0, ' [temp3 va] '));']);
        catch
            eval(['sol = struct2cell(solve(subs([grad_L ; us(ativas)], [' temp(1:end-1) '], [' temp2(1:end-1) '])==0, ' [temp3 va] '));']);
        end
        
        temp4 = sprintf('%f,',sol{end-length(variaveis)+1:end});
        temp5 = sprintf('%f,',sol{1:end-length(variaveis)});
        
        % Checa se todos us sao >= 0
        if ~all(cellfun(@isempty,sol)) && all([sol{1:end-length(variaveis)}] >= -1e-3)
            % Checa se todas as restricoes sao satisfeitas
            if eval(['all(subs(us, [' va ',' temp(1:end-1) ',' temp3(1:end-1) '], [' temp4(1:end-1) ',' temp2(1:end-1) ',' temp5(1:end-1) ']) <= 1e-3)'])
                % Esta solucao satisfaz as condições de KKT de 1a ordem
                eval(['Sols(:,k) = [' temp4(1:end-1) ',' temp2(1:end-1) ',' temp5(1:end-1) ']']);
                k = k + 1;
            end
        end
    elseif length(inativas) > 0 && length(ativas) == 0
        % Todas inativas
        temp = sprintf('mu%d,',inativas);
        temp2 = sprintf('%d,',zeros(1,length(inativas)));
        
        try
            eval(['sol = struct2cell(vpasolve(subs(grad_L, [' temp(1:end-1) '], [' temp2(1:end-1) '])==0, ' va '));']);
        catch
            eval(['sol = struct2cell(solve(subs(grad_L, [' temp(1:end-1) '], [' temp2(1:end-1) '])==0, ' va '));']);
        end
        
        temp4 = sprintf('%f,',sol{end-length(variaveis)+1:end});
        
        % Checa se todas as restricoes sao satisfeitas
        if ~all(cellfun(@isempty,sol)) && eval(['all(subs(us, [' va ',' temp(1:end-1) '], [' temp4(1:end-1) ',' temp2(1:end-1) ']) <= 1e-3)'])
            % Esta solucao satisfaz as condições de KKT de 1a ordem
            eval(['Sols(:,k) = [' temp4(1:end-1) ',' temp2(1:end-1) ']']);
            k = k + 1;
        end
    elseif length(inativas) == 0 && length(ativas) > 0
        % Todas ativas
        temp3 = sprintf('mu%d,',ativas);
        
        try
            eval(['sol = struct2cell(vpasolve([grad_L ; us], ' [temp3 va] '));']);
        catch
            eval(['sol = struct2cell(solve([grad_L ; us], ' [temp3 va] '));']);
        end
        
        temp4 = sprintf('%f,',sol{end-length(variaveis)+1:end});
        temp5 = sprintf('%f,',sol{1:end-length(variaveis)});
        
        % Checa se todos us sao >= 0
        if ~all(cellfun(@isempty,sol)) && all([sol{1:end-length(variaveis)}] >= -1e-3)
            % Checa se todas as restricoes sao satisfeitas
            if eval(['all(subs(us, [' va ',' temp3(1:end-1) '], [' temp4(1:end-1) ',' temp5(1:end-1) ']) <= 1e-3)'])
                % Esta solucao satisfaz as condições de KKT de 1a ordem
                eval(['Sols(:,k) = [' temp4(1:end-1) ',' temp5(1:end-1) ']']);
                k = k + 1;
            end
        end
    end
end

for i = 1:size(Sols,2)
    fprintf('\nSolucao %d:\n', i);
    
    fprintf(['x* = [%.6f %.6f %.6f]\n'...
        'mu* = [%.6f %.6f %.6f %.6f %.6f %.6f %.6f]\n'], ...
        Sols(:,i).');
    
    fprintf('Restricoes de desigualdade: [%.6f %.6f %.6f %.6f %.6f %.6f %.6f]\n', ...
        subs(us, [x1,x2,x3,mu1,mu2,mu3,mu4,mu5,mu6,mu7], Sols(:,i).'));
    
    fsol = subs(fo, [x1,x2,x3], Sols(1:3,i).');
    Lsol = subs(L, [x1,x2,x3,mu1,mu2,mu3,mu4,mu5,mu6,mu7], Sols(:,i).');
    fprintf('L(x*) = %.6f\n', Lsol);
    fprintf('f(x*) = %.6f\n', fsol);
    fprintf('p* - d* = %.6f\n\n', fsol - Lsol);
end