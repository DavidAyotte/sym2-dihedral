
print("read : bernoulli.gp")

isint(n) = {
	return(n == floor(n));
}
addhelp(isint, "isint(n) : vérifie si n est un entier ou non");

countdigit(n) = {
	if(!isint(n), error("n doit être un entier"));
	n = abs(n);
	i = 1;
	while(n > 10^i, i = i+1);
	return(i);
addhelp(countdigit, "countdigit(n) : calcul le nombre de chiffre d'un nombre");}

berndenom(n) = {
	if(!isint(n), error("n doit être un entier"));
	my(denom = 1);
	for(p=1,n+1, if(isprime(p) && n%(p-1) == 0, denom *= p));
	return(denom);
}
addhelp(berndenom, "denominateur du n-ieme nombre de bernoulli")

unsbernnum(n) = {
	if(n == 0, return(0));
	if(n == 1, return(1));
	if(n >= 3 && n%2 == 1, return(0));
	d = berndenom(n);
	prec = ceil(n*log(n) + log(d));
	default(realprecision, prec);
	K = (2*n!)/(2*Pi)^(n);
	M = ceil((K*d)^(1/(n-1)));
	z = prodeuler(p = 1, M, (1 - p^(-n))^(-1));
	return(ceil(d*K*z));
}
addhelp(unsbernnum, "Numerateur du n-ieme nombre de bernoulli")

