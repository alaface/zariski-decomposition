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
