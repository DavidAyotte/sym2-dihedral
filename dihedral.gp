
print("read : dihedral.gp");

mfispordinary(f, p) = {
	a = mfcoef(f, p);
	Na = norm(a);
	if(Na == 0, return(0));
	return(1);
}

dih(f, nbf, p, maxcoefs=150) = {
	forprime(q = 1, maxcoefs,
		if (kronecker(p, q) == -1,
			c = norm(mfcoef(f, q));
			if (c%p != 0,return(0), next(););
		,
			next();
		);
	);
	return(1);
};

exidih(S, p, maxcoefs=150) = {
	Base = mfeigenbasis(S);
	L = length(Base);
	nbf = nfinit(mffields(S)[1]);
	if (L>0,
		for(i = 1, L,
			f = Base[i];
			if (mfisdihedral(f, p, nbf, maxcoefs=maxcoefs), return(1),next())
		);
	,
		return(0)
	);
	return(0);
};

mfisdihedral(f, p, {D = 1}, {maxcoefs=150}) = {
	if(D == 1, D = p; chiD(q) = kronecker(q, D), chiD(q) = kronecker(D, q));
	forprime(q = 1, maxcoefs,
		if (chiD(q) == -1,
			cf = norm(mfcoef(f, q));
			if (cf%p != 0, return(0));
		);
	);
	return(1);
};
addhelp(mfisdihedral, "mfisdihedral(f, p, maxcoefs=150) : Verifie si les maxcoefs premier coefficients de la forme modulaire f satisfait la condition de congruence dihedrale pour un nombre premier p");

mfexistdihedral(S, p, maxcoefs=150) = {
	Base = mfeigenbasis(S);
	L = length(Base);
	if (L>0,
		for(i = 1, L,
			f = Base[i];
			if (mfisdihedral(f, p, maxcoefs=maxcoefs), return(1),next())
		);
	,
		return(0)
	);
	return(0);
};
addhelp(mfexistdihedral, "mfexistdihedral(S, p) : Determine s'il existe une eigenform dans l'espace S qui satisfait la congruence dihedral pour un nombre premier p");

testdihedral(lc = 1, pmax = 167, pmod = 3) = {
	forprimestep(p = pmod, pmax, 4,
		k = (p+1)/2;
		l = k + (lc)*(p-1);
		S = mfinit([1,l], 1);
		test = mfexistdihedral(S, p);
		if(test,
			printf("p = %s, k = %s, l = %s\nCongruence dihedrale trouvee !\n\n", p, k, l)
		,
			printf("p = %s, k = %s, l = %s\nAucune congruence dihedrale\n\n", p, k, l)
		);
	);
};
addhelp(testdihedral, "testdihedral(lc = 1, pmax = 167, pmod = 3) : Test la congruence dihedral pour plusieurs espaces");



