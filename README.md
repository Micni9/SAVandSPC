# Numerical Experiment
## Example 1

We consider the following equation:
$$
    \ut = \Delta u - h(u),x\in\Omega,t\in[0,T]
$$
Here $\Delta$ is the Laplace operator. $\Omega=[0,2\pi]$
Here the nonlinear function $h(u)=100\ f(u)$ with the following definition:
$$
    f(u) =\left\lbrace  
    \begin{array}{ll}
     2(u-1),    & u >1 \\
     (u^{3}-u), & u\in[-1, 1] \\
     2(u+1),    & u< -1
    \end{array}
    \right.
$$
The nonlinear function $f$ is modified from $f_0(u)=(u^3-u)$ to enforce the Lipschitz condition on $f$ with Lipschitz constant $2$, thus $h$ is also Lipschitz with the Lipschitz constant:
$$
    L=\max_{|u|\leq1}|h'(u)|=200
$$
This equation can be solved by both methods that we cover. We compare the stability performance of Stabilized Predictor-Corrector Method with the SAV\cite{ShenJie2019ANCo} and its variants for different $\frac{\Delta x}{\dt}= 2\pi$ by computing the energy functional $\varepsilon(u)=\int_\Omega\frac12\nrm{\nabla u}+25(u^2-1)^2dx$:\\
For SAV and its variants, $\EEEE(u) = \int_\Omega25(u^2-1)^2dx$
In this example, we choose the initial value $u(x,0) = \sin x,T=20$\\
