# Numerical Experiment
## Example 1

We consider the following equation:
$$ u_t = \Delta u - h(u),x\in\Omega,t\in[0,T] $$
Here $\Delta$ is the Laplace operator. $\Omega=[0,2\pi]$
Here the nonlinear function $h(u)=100\ f(u)$ with the following definition:
$$
\begin{equation*} 
f(u) =\left\lbrace  
    \begin{array}{ll}
     2(u-1),    & u >1 \\
     (u^{3}-u), & u\in[-1, 1] \\
     2(u+1),    & u< -1
    \end{array}
    \right.
\end{equation*}
$$
The nonlinear function $f$ is modified from $f_0(u)=(u^3-u)$ to enforce the Lipschitz condition on $f$ with Lipschitz constant $2$, thus $h$ is also Lipschitz with the Lipschitz constant:
$$
    \begin{equation*}L=\max_{|u|\leq1}|h'(u)|=200\end{equation*}
$$
This equation can be solved by both methods that we cover. We compare the stability performance of Stabilized Predictor-Corrector Method with the SAV and its variants for different $\frac{\Delta x}{\Delta t}= 2\pi$ by computing the energy functional $\varepsilon(u)=\int_\Omega\frac12\|{\nabla u}\|^2+25(u^2-1)^2dx$:\
For SAV and its variants, $\varepsilon_1(u) = \int_\Omega25(u^2-1)^2dx$
In this example, we choose the initial value $u(x,0) = \sin x,T=20$\
