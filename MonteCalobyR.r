a=2;
b=5;
g=function(x)
{
  exp(-x^2/2)/sqrt(2*pi);
}
f=function(y)
{
  (g(a+(b-a)*y)-c)/(d-c);
}
c=g(5);d=g(2);s_0=(b-a)*(d-c);
n=10^4;
x=runif(n);y=runif(n);
mu_n=sum(y<=f(x));
J=mu_n/n;
J_0=s_0*J+c*(b-a);
J_0
