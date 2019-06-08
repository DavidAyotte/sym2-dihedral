
print("read : dihedral.gp");

mfispordinary(f, p) = {
	a = mfcoef(f, p);
	Na = norm(a);
	if(Na == 0, return(0));
	return(1);
}
addhelp(mfispordinary, {"mfispordinary(f, p) : vérifie si la norme du p-ième coefficient de f est divisible par p.\n\n
Paramètres :\n
- f : une forme modulaire;\n 
- p : un premier rationnel."})

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
addhelp(mfisdihedral, {"mfisdihedral(f, nf, p, primeover, {D = 1}, {maxcoefs=150}) : Vérifie si une forme modulaire satisfait une congruence diédrale en p.\n\n
Paramètres :\n
- f : une forme modulaire;\n
- nf : un corps de nombre;\n
- p : un premier rationnel;\n
- primeover : un premier de nf au dessus de p;\n 
- D : un entier congru à 0 ou 1 mod 4. Si D = 1, alors la fonction vérifie si la forme est diédrale au sens de Dummigan-Heim. Si D>1, alors la fonction vérifie si la forme est diédrale
au sens de Brown-Ghate;\n 
- maxcoefs : seuil maximal de coeffients à vérifier.\n"});

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
addhelp(mfexistdihedral, {"mfexistdihedral(S, p, {maxcoefs=150}) : Détermine s'il existe une forme cuspidale propre normalisée qui satisfait une congruence diédrale en p au sens de Dummigan-Heim.\n\n 
Paramètres :\n
- S : un espace de formes cuspidales;\n
- p : un premier rationnel;\n 
- maxcoefs : seuil maximal de coeffients à vérifier.\n"})


/*
Vérification de la congruence diédrale au sens de Dummigan-Heim des espaces de formes cuspidales de poids k = (p+1)/2 pour p entre 3 et 167.
Attention : ce test est très long à exécuter.
*/
test = 0; \\Pour exéctuer le test, changer ce paramètre à 1.
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

