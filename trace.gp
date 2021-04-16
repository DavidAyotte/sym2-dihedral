
print("read : trace.gp");

read("misc.gp");

facteurck(k) = {
	return(2^(2*k-3)*(k/2 - 1)!/(k-1)!);
}

traceniv1(k) = {
	if (k%2 != 0, error("k doit etre pair !"));
	ck = facteurck(k);
	bernoullipol = bernpol(k-1);
	
	\\ Nombres de Bernoulli generalises
	bern4 = berngen(k, 4);
	bern3 = berngen(k, 3);
	
	return(uk*(bern4 + 2*bern3 + bernfrac(2*k-2)*(1+k/bernfrac(k))));
}
addhelp(traceniv1, "traceniv1(k) : compute the trace of the symmetric square L-function");

ordpfactorial(n, p) = {
	if(!isprime(p), error("p doit être un nombre premier"));
	m = log(n)/log(p);
	if(floor(m) == m, 
		return(1),
		return(sum(i = 1, floor(m), floor(n/p^i)))
	);
};
addhelp(ordpfactorial, "Calcul efficace de la valuation p-adique de n!");

ordptraceniv1(k, p) = {
	if(!isprime(p) || p==2, error("p doit être un nombre premier > 2"));
	if(k%2 != 0, error("k doit etre pair !"));
	bern4 = berngen(k, 4);
	bern3 = berngen(k, 3);
	B = bern4 + 2*bern3 + bernfrac(2*k-2)*(1+k/bernfrac(k));
	ordp = ordpfactorial(k/2 - 1, p) - ordpfactorial(k - 1, p) + valuation(B, p);
	return(ordp);
}
addhelp(ordptraceniv1, "p-valuation de la trace de la fonction L carrée symétrique")

bernsum(k) = {
	bern4 = berngen(k, 4);
	bern3 = berngen(k, 3);
	B = bern4 + 2*bern3 + bernfrac(2*k-2)*(1+k/bernfrac(k));
	return(B);
}
addhelp(bernsum, "Calcul du facteur comprenant une somme de nombres de Bernoulli dans la formule de la trace")

vpbernsum(k, p) = {
	return(valuation(bernsum(k), p));
}
addhelp(vpbernsum, "p-valuation de bernsum")

checkmp(p, m) = {
	if((2*m+1)%p == 0, return(1));
	return(0);
}
addhelp(checkmp, "vérification rapide")
