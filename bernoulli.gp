
isint(n) = {
	return(n == floor(n));
}
countdigit(n) = {
	if(!isint(n), error("n doit être un entier"));
	n = abs(n);
	i = 1;
	while(n > 10^i, i = i+1);
	return(i);
}
berndenom(n) = {
	if(!isint(n), error("n doit être un entier"));
	my(denom = 1);
	for(p=1,n+1, if(isprime(p) && n%(p-1) == 0, denom *= p));
	return(denom);
}
bernnum1(n) = {
	if(!isint(n), error("n doit être un entier"));
	my(denom = berndenom(n));
	return(bernfrac(n)*denom);
}

\\forstep(n = 2, 1000, 2,{
\\	printf("n = %s,  nb = %s\n", n, countdigit(bernnum1(n)));
\\});

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

/* p = 167;
k = (p+1)/2;
l = k+300*(p-1);
denom = 1;
for(p = 1, l,{
	if(isprime(p) && l%(p-1)==0, denom = p*denom);
});
\\print(denom)
B = bernfrac(l);
num = -floor(denom*B);
\\if(floor(num) == num, print("ok"));
print(num%10000000);
\\i = 1;
\\while(1, {
\\	if(num < 10^i,print(num < 10^i); print(10^i); print(i); break);
\\	i = i+1;
\\}) */

