
print("read : dihedral.gp");
addhelp(dihedral, {"
This gp file contains two functions about dihedral congruences of modular form. 

Let p be a prime and k an even integer. Let f be a primitive Hecke eigenform with coefficients a_n(f). Let p be a ratonal prime and P | p a prime of the Hecke field of f.

We say that f satisfies a dihedral congruence at p in the sense of Dummigan-Heim if:
	
(DH) P | a_q(f) for every prime q s. t. (q/p) = -1, where (./p) is the Kronecker symbol.
	
Similarly, we say that f satisfies a dihedral congruence at p in the sense of Brown-Ghate if
	
(BG) P | (D/q)*a_q(f) for every prime q s.t. (D/q) = -1, where D is an integer congruent to 0 or 1 mod 4.

List of functions :
- mfisdihedral : Return 1 if f satisfies a dihedral congruence in the sense of DH or BG at p and 0 if not;
- mfexistdihedral : Return 1 if there exists a modular form in a cuspidal subspace that satisfies a dihedral congruence in the sense of Dummigan-Heim and 0 if not.

Type ?? <function name> for more detail 

References
- Ayotte, D, Relations entre le nombre de classes et les formes modulaires, master thesis (2019)
- Dummigan, Neil; Heim, Bernhard, Symmetric square L-values and dihedral congruences for cusp forms. J. Number Theory 130 (2010)
- Brown, Alexander F.; Ghate, Eknath P. Dihedral congruence primes and class fields of real quadratic fields. J. Number Theory 95 (2002)
"})

mfisdihedral(f, nf, p, primeover, {D = 1}, {maxcoefs=150}) = {
	if(D == 1, D = p; chiD(q) = kronecker(q, D), chiD(q) = kronecker(D, q));
	forprime(q = 1, maxcoefs,
		if (chiD(q) == -1,
			cf = mfcoef(f, q);
			if (nfeltreduce(nf, cf, primeover) != 0, return(0));
		);
	);
	return(1);
};
addhelp(mfisdihedral, {"mfisdihedral(f, nf, p, primeover, {D = 1}, {maxcoefs=150}) : Return 1 if f satisfies a dihedral congruence at p and 0 if not.
Parameters
- f : a modular form;
- nf : a number field;
- p : a rational prime;
- primeover : a prime of nf over p;
- D : an integer congruent to 0 or 1 mod 4. If D = 1, then the function verifies whether f is dihedral in the sense of Dummigan-Heim. If D>1, then the function verifes whether f is dihedral in the sense of Brown-Ghate;
- maxcoefs : the maximal index to verify the congruence"});

mfexistdihedral(S, p, {maxcoefs=150}) = {
	nfs = mffields(S);
	B = mfeigenbasis(S);
	if(#B == 0, return(0));
	for(i = 1, #B,
		f = B[i];
		fld = nfinit(nfs[i]);
		P = idealprimedec(fld, p);
		for(j = 1, #P,
			primeover = idealhnf(fld, P[j]);
			if (mfisdihedral(f, fld, p, primeover, D=1, maxcoefs=maxcoefs), return(1));
		);
	);
	return(0);
};
addhelp(mfexistdihedral, {"mfexistdihedral(S, p, {maxcoefs=150}) : Return 1 if there exists a modular form in S that satisfies a dihedral congruence in the sense of Dummigan-Heim.
Parameters
- S : a cuspidal subspace;
- p : a rational prime; 
- maxcoefs : the maximal index to verify the congruence"})


/*
Check dihedral congruence in the sense of Dummigan-Heim for the cuspidal subspaces of weight k = (p+1)/2 for p between 3 and 167
Caution : this test is very long to run !
*/
test = 0; \\To execute this test, change this paramter to 1.
if(test,{
	forprimestep(p = 3, 167, 4,
		k = (p+1)/2;
		S = mfinit([1,k], 1);
		if(mfexistdihedral(S, p),
			printf("p = %s, k = %s\nDihedral congruence found !\n\n", p, k);
		,
			printf("p = %s, k = %s\nNo dihedral congruence\n\n", p, k);
		);
	);
});

