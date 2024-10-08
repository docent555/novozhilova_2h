function res = gmp(m,n,Rb,Rr,estx)

syms x
f = diff(besselj(m,x),x);
% df=inline(f);
df=matlabFunction(f);
nu = fzero(df, estx);
z=0:0.01:100;
plot(z,df(z))
hold on
plot(nu,df(nu),'ro')

res = (besselj(m-n, nu*Rb/Rr)/ besselj(m, nu))^2/(nu^2 - m^2);

end

