function [eta1, eta2] = eff(eta1, eta2, f, Dtr, Ne, Nt, Zex, u, n) %#codegen

fp = complex(zeros(2,1));

th0 = 2.0D0*pi*(0:Ne-1)/Ne;
p0 = exp(1i*th0).';

opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

for i=1:Nt
    fp(1) = f(i,1)*exp(1i*f(i,2));
    fp(2) = f(i,3)*exp(1i*f(i,4));

    p1 = oscill_ode(fp(1), Ne, [0 Zex], Dtr(1), p0, u, n, opts);
    p2 = oscill_ode(fp(2), Ne, [0 Zex], Dtr(2), p0, u, n, opts);

    eta1(i) = 1 - mean(abs(p1(end,:)).^2);
    eta2(i) = 1 - mean(abs(p2(end,:)).^2);
end
end

