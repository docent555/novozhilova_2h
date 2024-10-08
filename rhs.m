function s = rhs(z, p, delta, a, u, n)

s = 1i*(a*u(z)*conj(p).^(n - 1) - (delta + abs(p).^2 - 1).*p);

end
