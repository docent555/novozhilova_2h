function r = xi(p, n, u, zax)

m = conj(u(zax)).'.*mean(p.^n, 2);
r = trapz(zax, m);

end

