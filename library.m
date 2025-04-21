///////////////////////////////////////////////////////////////////////////
// Support functions for Zariski decomposition over a rational lattice
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// IntersectionNumber(D, C, M)
// Input: 
//   - D: vector representing a divisor
//   - C: vector representing a curve
//   - M: symmetric intersection matrix
// Output:
//   - The intersection number D·C computed as Transpose(D)*M*C
///////////////////////////////////////////////////////////////////////////
IntersectionNumber := function(D, C, M)
    return (D*M*C)[1];
end function;

///////////////////////////////////////////////////////////////////////////
// NegativeCurves(D, Curves, M)
// Input:
//   - D: divisor vector
//   - Curves: list of curve vectors
//   - M: intersection matrix
// Output:
//   - List of curves C in Curves such that D·C < 0
///////////////////////////////////////////////////////////////////////////
NegativeCurves := function(D, Curves, M)
    return [ C : C in Curves | IntersectionNumber(D, C, M) lt 0 ];
end function;

///////////////////////////////////////////////////////////////////////////
// SolveNegativePart(D, NegC, M)
// Input:
//   - D: divisor vector
//   - NegC: list of curves with D·C < 0
//   - M: intersection matrix
// Output:
//   - List of rational coefficients b_i such that 
//     N = Σ b_i C_i satisfies N·C_i = D·C_i
///////////////////////////////////////////////////////////////////////////
SolveNegativePart := function(D, NegC, M)
    q := #NegC;
    if q eq 0 then return []; end if;

    A := ZeroMatrix(Rationals(), q, q);
    b := ZeroMatrix(Rationals(), q, 1);
    for i in [1..q] do
        for j in [1..q] do
            A[i,j] := IntersectionNumber(NegC[i], NegC[j], M);
        end for;
        b[i,1] := IntersectionNumber(D, NegC[i], M);
    end for;
    sol := Solution(A, b);
    return [ sol[i,1] : i in [1..q] ];
end function;

///////////////////////////////////////////////////////////////////////////
// IsNef(D, Curves, M)
// Input:
//   - D: divisor vector
//   - Curves: list of curve vectors
//   - M: intersection matrix
// Output:
//   - true if D·C ≥ 0 for all C in Curves (i.e., D is nef), false otherwise
///////////////////////////////////////////////////////////////////////////
IsNef := function(D, Curves, M)
    return &and[ IntersectionNumber(D, C, M) ge 0 : C in Curves ];
end function;

///////////////////////////////////////////////////////////////////////////
// ZariskiDecomposition(D, Curves, M)
// Input:
//   - D: pseudo-effective divisor vector
//   - Curves: list of irreducible curve vectors
//   - M: intersection matrix
// Output:
//   - A pair <P, N> such that D = P + N is the Zariski decomposition,
//     where P is nef and N is effective with negative definite support
///////////////////////////////////////////////////////////////////////////
ZariskiDecomposition := function(D, Curves, M)
    V := Parent(D);
    P := D;
    Ntot := V!0;

    repeat
        Neg := NegativeCurves(P, Curves, M);
        if #Neg eq 0 then
            return P, Ntot;
        end if;

        coeff := SolveNegativePart(P, Neg, M);
        Nstep := V!0;
        for i in [1..#Neg] do
            Nstep +:= coeff[i]*Neg[i];
        end for;

        P    -:= Nstep;
        Ntot +:= Nstep;
    until IsNef(P, Curves, M);

    return P, Ntot;
end function;

// ZariskiDecompositionCox
// Input:
//   - wD: an element of a rational vector space, representing the class of a divisor D
//   - W: a list [w1, ..., wr] of vectors in the same space, representing the degrees of the generators f_i of Cox(X)
//
// Output:
//   - wP: the class of the positive part P of the Zariski decomposition of D
//   - wN: the class of the negative part N = ∑ μ_i w_i
//   - mu: the list of rational numbers [μ_1, ..., μ_r] such that D = P + ∑ μ_i D_i
//
// The algorithm computes μ_i by finding the minimal non-negative x such that 
// wD = x * w_i + y * R for some extremal ray R of the cone τ_i = cone(wD, -w_i) ∩ cone(w_j : j ≠ i)

ZariskiDecompositionCox := function(wD, W)
    r := #W;
    mu := [Rationals()!0 : i in [1..r]];

    for i in [1..r] do
        cone1 := Cone([wD, -W[i]]);
        cone2 := Cone([W[j] : j in [1..r] | j ne i]);
        tau := cone1 meet cone2;
        rays := Rays(tau);

        candidates := [];

        for R in rays do
            M := Matrix(Rationals(), [Eltseq(W[i]), Eltseq(R)]);
            success, sol := IsConsistent(M, Vector(wD));
            if success then
                x := sol[1];
                if x ge 0 then
                    Append(~candidates, x);
                end if;
            end if;
        end for;

        mu[i] := Minimum(candidates);
    end for;

    wN := &+[mu[i]*W[i] : i in [1..r]];
    wP := wD - wN;

    return wP, wN, mu;
end function;
