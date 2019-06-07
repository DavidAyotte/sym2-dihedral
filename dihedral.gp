
print("read : dihedral.gp");

mfispordinary(f, p) = {
	a = mfcoef(f, p);
	Na = norm(a);
	if(Na == 0, return(0));
	return(1);
}

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
addhelp(mfisdihedral, "mfisdihedral(f, p, {D=1}, {maxcoefs=150}) : Verifie si les maxcoefs premier coefficients de la forme modulaire f satisfait la condition de congruence dihedrale pour un nombre premier p");

mfexistdihedral(S, p, maxcoefs=150) = {
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

test = 0;
if(test,{
	forprimestep(p = 3, 167, 4,
		k = (p+1)/2;
		S = mfinit([1,k], 1);
		if(mfexistdihedral(S, p),
			printf("p = %s, k = %s\nCongruence dihedrale trouvee !\n\n", p, k);
		,
			printf("p = %s, k = %s\nAucune congruence dihedrale\n\n", p, k);
		);
	);
});

