# Numerical Methods for MACM 316 Students

Below is a list of the available numerical methods.

## Linear Systems of Equations: Iterative Methods

- Jacobi method (`jacobi.m`), based on [BF, Algorithm 7.1].
- Gauss-Seidel method (`gauss_seidel.m`), based on [BF, Algorithm 7.2].
- Successive Over-Relaxation (`sor.m`), based on [BF, Algorithm 7.3].

## Nonlinear Equations

- Bisection method (`bisection.m`) based on [BF, Algorithm 2.1].
- Newton's method (a.k.a. Newton-Raphson's method, `newton_raphson.m`) based on [BF, Algorithm 2.3].
- Secant method (`secant.m`) based on [BF, Algorithm 2.4].

### Demos
- `demo_bisection.m` Solves a nonlinear equation using bisection method.
- `demo_bisection_newton.m` Solves three different equations using bisection, Newton's, and secant method.
- `demo_multiple_roots.m` Illustrates the performance of Newton's method when the roots have multiplicity > 1.

## References

[BF] Burden, R. L., Faires, D. J., & Burden, A. M. (2015). Numerical Analysis, 10th Edition. Cengage Learning.
