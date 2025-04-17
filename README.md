# Zariski Decomposition in Magma

This repository provides a Magma implementation of the **Zariski decomposition** of a pseudo-effective divisor on an algebraic surface, based on its intersection matrix and a list of irreducible curve classes.

## Contents

- `library.m` — Core library implementing the decomposition algorithm.
- `example_dim3.m` — Example in dimension 3 with a tridiagonal intersection matrix.

## Description

Given a surface $X$, a pseudo-effective divisor $D$, and a finite set of irreducible curves $C_1, \dots, C_n$, the **Zariski decomposition** expresses $D$ uniquely as:

$$
D = P + N
$$

where:
- $P$ is a **nef** divisor (i.e. $P \cdot C_i \ge 0$ for all $i$),
- $N = \sum b_i C_i$ is an **effective** divisor with support on a negative definite intersection matrix,
- $P \cdot N = 0$.

## Mathematical Reference

This algorithm is based on:

**Oscar Zariski**,  
*The Theorem of Riemann-Roch for High Multiples of an Effective Divisor on an Algebraic Surface*,  
Annals of Mathematics, Vol. 76, No. 3 (Nov., 1962), pp. 560–615.  
See **Theorem 7.7** for the step-by-step construction.

## How it Works

The algorithm iteratively removes the negative part of $D$ as follows:

1. **Negative curves**: Identify all curves $C_i$ such that $D \cdot C_i < 0$.
2. **Linear system**: Solve the system  
   $$
   (C_i \cdot C_j)_{i,j} \cdot \mathbf{b} = (D \cdot C_i)_i
   $$
where the $b_j$ are the coefficients of the negative part $N = \sum b_j C_j$.
3. **Update**: Replace $D := D - N$, and repeat.
4. **Termination**: The process stops when the updated divisor is nef. Since the Néron–Severi group has finite rank, this happens in finitely many steps.

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
