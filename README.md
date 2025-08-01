# KKT-Conditions-Tester-1st-and-2nd-Order-Matlab
MATLAB implementation for verifying 1st and 2nd order KKT optimality conditions. This tool automates the combinatorial check of all possible active and inactive constraint sets for first-order conditions and evaluates the Hessian of the Lagrangian for second-order conditions in constrained optimization problems.

---

### Condições de Karush-Kuhn-Tucker (KKT)
As condições de Karush-Kuhn-Tucker (KKT) são um conjunto de critérios necessários para que uma solução em um problema de otimização não-linear seja ótima. Elas generalizam o método dos multiplicadores de Lagrange para incluir restrições de desigualdade. Para um ponto ser considerado um mínimo local (sob certas condições de regularidade), ele deve satisfazer as condições KKT de primeira ordem. As condições de segunda ordem são então utilizadas para verificar se esse ponto é de fato um mínimo (ou máximo).

#### Exemplo Analisado

O seguinte problema de programação quadrática foi usado de exemplo para os scripts nesse repositório:

$$
\begin{aligned}
\text{Minimizar } f(x) = & 2x_1^2 + 20x_2^2 + 43x_3^2 + 12x_1x_2 - 16x_1x_3 - 56x_2x_3 \\
& + 8x_1 + 20x_2 + 6x_3 \\
\text{Sujeito a:} \\
& 3x_1 + 2x_2 + 5x_3 \le 35 \\
& x_1 + 2x_2 + 3x_3 \ge 5 \\
& -x_1 + 2x_2 - 5x_3 \le 3 \\
& 5x_1 - 3x_2 + 2x_3 \le 30 \\
& x_1, x_2, x_3 \ge 0
\end{aligned}
$$

---

### 1. Condições KKT de Primeira Ordem
As condições de primeira ordem estabelecem uma relação entre o gradiente da função objetivo e os gradientes das restrições ativas no ponto ótimo. Elas incluem:
1.  **Estacionariedade:** O gradiente do Lagrangiano deve ser nulo.
2.  **Viabilidade Primal:** Todas as restrições devem ser satisfeitas.
3.  **Viabilidade Dual:** Os multiplicadores de Lagrange associados às restrições de desigualdade devem ser não-negativos ($\mu_i \ge 0$).
4.  **Folga Complementar:** O produto de cada multiplicador de Lagrange pela sua respectiva restrição de desigualdade deve ser zero ($\mu_i g_i(x) = 0$).

Para escrever as condições necessárias de KKT, temos que derivar o Lagrangiano dado por:

$$
\begin{align}
&\mathcal{L} = \notag\\
&(2 \cdot x_1^2 + 20 \cdot x_2^2 + 43 \cdot x_3^2 + 12 \cdot x_1 \cdot x_2 - 16 \cdot x_1 \cdot x_3 - 56 \cdot x_2 \cdot x_3 + 8 \cdot x_1 + 20 \cdot x_2 + 6 \cdot x_3) \notag\\
&+ \mu_1 (3 \cdot x_1 + 2 \cdot x_2 + 5 \cdot x_3 - 35) \notag\\
&+ \mu_2 (-x_1 - 2 \cdot x_2 - 3 \cdot x_3 + 5) \notag\\
&+ \mu_3 (-x_1 + 2 \cdot x_2 - 5 \cdot x_3 - 3) \notag\\
&+ \mu_4 (5 \cdot x_1 - 3 \cdot x_2 + 2 \cdot x_3 - 30) \notag\\
&+ \mu_5 (-x_1) \notag\\
&+ \mu_6 (-x_2) \notag\\
&+ \mu_7 (-x_3)
\end{align}
$$

Assim, as condições necessárias de KKT de primeira ordem são:

$$
\begin{align}
&4 \cdot x_1 + 12 \cdot x_2 - 16 \cdot x_3 + 8 + 3 \cdot \mu_1 - \mu_2 - \mu_3 + 5 \cdot \mu_4 - \mu_5 = 0 \notag\\
&12 \cdot x_1 + 40 \cdot x_2 - 56 \cdot x_3 + 20 + 2 \cdot \mu_1 - 2 \cdot \mu_2 + 2 \cdot \mu_3 - 3 \cdot \mu_4 - \mu_6 = 0 \notag\\
&- 16 \cdot x_1 - 56 \cdot x_2 + 86 \cdot x_3 + 6 +  5 \cdot \mu_1 - 3 \cdot \mu_2 - 5 \cdot \mu_3 + 2 \cdot \mu_4 - \mu_7 = 0 \notag\\
&\mu_1 \geq 0, \mu_2 \geq 0, \mu_3 \geq 0, \mu_4 \geq 0, \mu_5 \geq 0, \mu_6 \geq 0, \mu_7 \geq 0 \notag\\
&\mu_1 (3 \cdot x_1 + 2 \cdot x_2 + 5 \cdot x_3 - 35) = 0  \notag\\
&\mu_2 (-x_1 - 2 \cdot x_2 - 3 \cdot x_3 + 5) = 0  \notag\\
&\mu_3 (-x_1 + 2 \cdot x_2 - 5 \cdot x_3 - 3) = 0  \notag\\
&\mu_4 (5 \cdot x_1 - 3 \cdot x_2 + 2 \cdot x_3 - 30) = 0  \notag\\
&\mu_5 (-x_1) = 0  \notag\\
&\mu_6 (-x_2) = 0  \notag\\
&\mu_7 (-x_3) = 0 
\end{align}
$$

Dessa forma, com 7 restrições de desigualdade teremos $2^{7} = 128$ possíveis combinações de restrições ativas e inativas. Uma pergunta: **você** gostaria de avaliar uma por uma manualmente? Existem casos em que o número de combinações é extrapolado para a casa dos milhares, inviabilizando qualquer operação manual.

#### Implementação: `first_order_KKT.m`
Este script **automatiza** completamente a verificação das condições KKT de primeira ordem. Em vez de testar manualmente todas as combinações, o código gera e resolve sistematicamente os sistemas de equações para cada caso, identificando todos os pontos que satisfazem as condições necessárias de otimalidade.

Neste caso, a solução fornecida pelo script é única:

$x^*$ = [0.285714 1.017857 0.892857]

$\mu^*$ = [0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 7.071429]

---

### 2. Condições KKT de Segunda Ordem
Uma vez que uma solução é encontrada pelas condições necessárias de KKT de primeira ordem, as condições de segunda ordem são usadas para classificá-lo. A condição suficiente de segunda ordem para um mínimo/máximo global exige que a matriz Hessiana do Lagrangiano seja definida positiva/negativa no subespaço do plano tangente $M$, contudo, se a Hessiana for semi-definida positiva/negativa, o ponto é classificado como ao menos um mínimo/máximo local. Para o exemplo, a Hessiana do Lagrangiano é:

$$
\mathcal{H}(L) = \begin{bmatrix}
4 & 12 & -16 \\
12 & 40 & -56 \\
-16 & -56 & 86
\end{bmatrix}
$$

A análise requer a projeção desta matriz no subespaço definido pelos gradientes das restrições ativas e a verificação dos autovalores da matriz resultante.

#### Implementação: `second_order_KKT.m`
Este script implementa a verificação da condição de segunda ordem. Ele testa a definição da Hessiana do Lagrangiano no subespaço do plano tangente $M$. O código realiza essa projeção e a análise de autovalores de forma automática, determinando a natureza de um ponto candidato a ótimo.

No exemplo, como é necessário que $L$ seja positivo definido no subespaço $M$ para que a solução seja mínima global, e no subespaço $M$ temos $y_1 = 0$, $y_2 = 0$ e $y_3 = 0$, $y^T L y = 0$, a solução não é minima global. Se não tivermos essas condições citadas, podemos reescrever a expressão:

$$
y^T \cdot L \cdot y = y_1\cdot (4\cdot y_1 + 12\cdot y_2 - 16\cdot y_3) + y_2\cdot (12\cdot y_1 + 40\cdot y_2 - 56\cdot y_3) - y_3\cdot (16\cdot y_1 + 56\cdot y_2 - 86\cdot y_3)
$$

$$
y^T \cdot L \cdot y = 4\cdot y_1^2 + 24\cdot y_1\cdot y_2 - 32\cdot y_1\cdot y_3 + 40\cdot y_2^2 - 112\cdot y_2\cdot y_3 + 86\cdot y_3^2
$$

Substituindo $M_2 = y_1 = - 2\cdot y_2 - 3\cdot y_3$:

$$
y^T \cdot L \cdot y = 8\cdot y_2^2 - 72\cdot y_2\cdot y_3 + 218\cdot y_3^2
$$

$$
y^T \cdot L \cdot y = \begin{bmatrix}
y_2 & y_3
\end{bmatrix} \cdot \begin{bmatrix}
8 & -36 \\
-36 & 218
\end{bmatrix} \cdot \begin{bmatrix}
y_2 \\
y_3
\end{bmatrix}
$$

Os autovalores da matriz são:

$$
\text{eig}\left( \begin{bmatrix}
8 & -36 \\
-36 & 218
\end{bmatrix} \right) = \begin{bmatrix}
2 \\
224
\end{bmatrix}
$$

Portanto, para $y_2 \neq 0$ e $y_3 \neq 0$, $y^T \cdot L \cdot y > 0$, o ponto seria ao menos mínimo local desconsiderando as condições citadas anteriormente. Ainda, o script `second_order_KKT.m` nos fornece a informação que a Hessiana de $L$ é semi-definida positiva em $M$, confirmando a hipótese.

---

## Referências

[1] [Bazaraa, 1993] M.S. Bazaraa, H.D. Sherali, C.M. Shetty, Nonlinear Programming, 2nd Ed., John Wiley, 1993.

[2] [Bertsekas, 1999] Dimitri P. Bertsekas. Nonlinear Programming. Athena Scientific. 2nd Ed., 1999.

[3] [Antoniou, 2007] Antoniou, A. and Lu, W.S. Practical Optimization: Algorithms and Engineering Applications, Springer, 2007

[4] [Vanderbei, 2012] Robert J. Vanderbei. Linear Programming: Foundations and Extensions, Fourth Ed., 2012, Springer

[5] [Luenberger 2008] D.G. Luenberger, Y. Ye, Linear and Nonlinear Programming, Third Ed., Addison Wesley, 2008, Springer
