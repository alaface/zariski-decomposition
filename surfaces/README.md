# Zariski Decomposition on surfaces

This repository provides a Magma implementation of the **Zariski decomposition** of an effective divisor on an algebraic surface, based on the intersection matrix of the irreducible curve in the support of $D$.

## Contents

- `library.m` — Core library implementing the decomposition algorithm.

## Description

Given a surface $X$, aneffective divisor $D$, and a finite set of irreducible curves $C_1, \dots, C_n$ which form the support of $D$, the **Zariski decomposition** expresses $D$ uniquely as:

$D = P + N$

where:
- $P$ is a **nef** divisor (i.e. $P \cdot C_i \ge 0$ for all $i$),
- $N = \sum b_i C_i$ is an **effective** divisor with support on a negative definite intersection matrix,
- $P \cdot N = 0$.

## Reference

This algorithm is based on **Theorem 7.7** of the following paper:

**Oscar Zariski**,  
*The Theorem of Riemann-Roch for High Multiples of an Effective Divisor on an Algebraic Surface*,  
Annals of Mathematics, Vol. 76, No. 3 (Nov., 1962), pp. 560–615.  

## How it Works

The algorithm computes the Zariski decomposition of a pseudo-effective divisor $D$ by iteratively extracting its negative part. It works as follows:

1. **Negative curves**: Identify all irreducible curves $C_i$ in the support of $D$ such that $D \cdot C_i < 0$.

2. **Construct the negative part**: Let $C_1, \dots, C_r$ be the curves found in step 1.  
   Look for a divisor of the form $N = \sum_{j=1}^r b_j C_j$ such that:

   $$N \cdot C_i = D \cdot C_i \quad \text{for all } i = 1, \dots, r$$

   The matrix $(C_i \cdot C_j)$ is negative definite, so the system has a unique solution in $\mathbb{Q}^r_{\geq 0}$.

3. **Update**: Set $D := D - N$. The new divisor satisfies $(D \cdot C_i = 0)$ for all curves in the support of $N$.

4. **Termination**: Repeat the process with the updated $D$. The algorithm stops when $D$ becomes nef, i.e., when $D \cdot C_i \ge 0$ for all curves.  
      
## How to Use

1. Load the library:

```magma
load "library.m";

## Example: Dimension 3

You can run the following example after loading the library:

```magma
load "library.m";

R := Rationals();
V := VectorSpace(R, 3);

// Negative definite intersection matrix
M := Matrix(R, 3, 3, [ -2, 1, 0,
                    1,-2, 1,
                    0, 1,-2 ]);

// Curve basis
C1 := V![1,0,0];
C2 := V![0,1,0];
C3 := V![0,0,1];
Curves := [C1, C2, C3];

// Divisor with non-negative coefficients
D := V![2,1,0];

// Compute Zariski decomposition
P, N := ZariskiDecomposition(D, Curves, M);

print "D =", D;
print "P =", P;
print "N =", N;

// Check intersection with each curve
for i in [1..#Curves] do
 printf "P · C%o = %o\n", i, IntersectionNumber(P, Curves[i], M);
end for;
