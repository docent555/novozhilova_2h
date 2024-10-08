function p = oscill_ode(f, Ne, ZAxis, Delta, p0, u, n, opts) %#codegen

% persistent i;
% if isempty(i)
%     i = 1;
% end

p0v = [real(p0); imag(p0)];
reidx = 1:Ne;
imidx = Ne+1:2*Ne;

% if i == 10000
%     for j=1:Nz
%         ss(j) =  S1(ZAxis(j));
%     end
%     figure;
%     plot(ZAxis, imag(ss), ZAxis, imag(field))
%     pause
% end

[~, pv] = ode45(@(z, pv) rhsv(z, pv, Delta, f, u, n, reidx, imidx) , ZAxis , p0v, opts);
p = pv(:,reidx) + 1i*pv(:,imidx);

% ip=[ZAxis' imag(p)];
% save test.dat ip -ascii
% pause

% i = i + 1;
end

