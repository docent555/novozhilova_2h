function zv = rhsv(z, pv, delta, a, u, n, reidx, imidx)

p = pv(reidx) + 1i*pv(imidx);

z = rhs(z, p, delta, a, u, n);

zv = [real(z); imag(z)];

end