# Performance Analysis of the Jacobi Iterative Method

This project investigates the **Jacobi iterative method** for solving a sequence of linear systems:

$$
A_n x = b_n
$$

where $A_n$ is a large, sparse, symmetric, positive definite matrix.  
The objective is to measure how the **number of iterations required for convergence** scales with the matrix size $n$ for two different matrix configurations and explain this scaling mathematically.

---

## System Generation

A specific family of **tridiagonal matrices** is used, originating from the discretization of a 1D differential equation.  

- **Matrix entries**:  
  - Main diagonal: $c/h^2$, where $h = 1/(n+1)$  
  - Off-diagonals: $-1/h^2$  

- **Right-hand side vector**: $b_n$ is a vector of ones.  
- **Exact solution**: Computed using a direct method to serve as the "ground truth" for error measurement.

---

## Phase A — Convergence Analysis for Parameter $c=2$

This phase examines the Jacobi method performance for **$c=2$**.

### Algorithm Implementation

For the system $Ax = b$, the Jacobi iterative update for each element $x_i$ is:

$$
x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j \ne i} a_{ij} x_j^{(k)} \right)
$$

- For tridiagonal matrices, the sum reduces to only **two non-zero terms**.

### Experimental Procedure

- For increasing matrix sizes $n$:  
  - Generate $A_n$ and $b_n$ with $c=2$.  
  - Start Jacobi iteration from $x^{(0)} = 0$.  
  - Continue until the **relative error in the $\ell_\infty$ norm**:

$$
\frac{\|x^{(k)} - x\|_\infty}{\|x\|_\infty} < \text{TOL} = 10^{-3}
$$

- Record the **number of iterations** required for convergence.

### Complexity Analysis

- Plot **iterations vs $n$** on a **log-log scale**.  
- Slope of the line, $m$, reveals scaling:

$$
\text{Iterations} \propto n^m
$$

**Expected Result:** $m = 2$, since the number of iterations is proportional to the **matrix condition number**, which scales as $O(n^2)$.

---

## Phase B — Convergence Analysis for Parameter $c=4$

This phase repeats the experiment with **$c=4$**.

### Experimental Procedure

- Follow the same steps as in Phase A, but generate $A_n$ with $c=4$.  
- Record the number of iterations for each $n$.

### Complexity Analysis and Comparison

- Plot iterations vs $n$ on a **log-log scale**, along with Phase A results.  
- The slope remains $m = 2$.  
- Changing $c$ affects the **constant factor** in the number of iterations but **does not change the quadratic scaling** with $n$.

---

## What This Investigation Demonstrates

- **Hands-on implementation of Jacobi method:** Moves from textbook formula to working procedure.  
- **Numerical convergence analysis:** Measures **iterations to achieve desired accuracy**.  
- **Empirical scalability:** Shows how iteration count grows with problem size; confirms $O(n^2)$ scaling.  
- **Order vs constant factor:** Changing $c$ affects absolute iterations but not the growth rate.  
- **Log-log plots for power laws:** Visualizes and quantifies scaling, a key technique in scientific computing.

---
