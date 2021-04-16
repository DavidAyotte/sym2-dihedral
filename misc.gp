
print("read : misc.gp")

addhelp(misc, {"
This gp file contains multiple miscellaneous functions about Bernoulli numbers and modular forms.

List of functions :
- isint(n) : determine whether n is an integer of not;
- countdigitint(n) : count the number of digits of the integer n;
- berngen(k, D) : k-th generalised Bernoulli number for the Kronecker character (D/.);
- berndenom(n) : return the unsigned denominator of the n-th Bernoulli number;
- unsbernnum(n) : return the unsigned numerator of the n-th Bernoulli number;
- mfispordinary(f, p) : verify wheher the norm of the p-th coefficient of f is divisible by p.
"})

isint(n) = {
	return(n == floor(n));
}
addhelp(isint, "isint(n) : determine if n is an integer of not.");

countdigitint(n) = {
	if(!isint(n), error("n doit être un entier"));
	n = abs(n);
	i = 1;
	while(n > 10^i, i = i+1);
	return(i);
}
addhelp(countdigitint, "countdigitint(n) : count the number of digits of the integer n");

berngen(k, D) = {
	bernoullipol = bernpol(k-1);
	return D^(k-2)*sum(i = 1, D, kronecker(-D, i)*subst(bernoullipol, x, i/D));
}
addhelp(berngen, {"berngen(k, D) : k-th generalised Bernoulli number for the Kronecker character (D/.)"});

berndenom(n) = {
	if(n == 0 || (n >= 3 && n%2 == 1), return(1));
	if(n == 1, return(2));
	if(!isint(n), error("n doit être un entier"));
	my(denom = 1);
	for(p=1,n+1, if(isprime(p) && n%(p-1) == 0, denom *= p));
	return(denom);
}
addhelp(berndenom, "berndenom(n) : return the unsigned denominator of the n-th Bernoulli number");

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
addhelp(unsbernnum, "unsbernnum(n) : return the unsigned numerator of the n-th Bernoulli number");

mfispordinary(f, p) = {
	a = mfcoef(f, p);
	Na = norm(a);
	if(Na == 0, return(0));
	return(1);
}
addhelp(mfispordinary, {"mfispordinary(f, p) : verify whether the norm of the p-th coefficient of f is divisible by p.
Parameters :
- f : a modular form;
- p : a rationnal prime."})

