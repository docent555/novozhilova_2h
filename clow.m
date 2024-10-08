function [cl1, cl2] = clow(cl1, cl2, eta1, eta2, f, Th, R, n, I) %#codegen

lhs1 = 2*n*f(:,1).^2 - 4*n*R(1)*f(:,5).*f(:,1).*cos(Th(1) - f(:,6) + f(:,2));
rhs1 = I(1)*eta1;
cl1(:) = lhs1 - rhs1;

lhs2 = 2*n*f(:,3).^2 - 4*n*R(2)*f(:,5).*f(:,3).*cos(Th(2) - f(:,6) + f(:,4));
rhs2 = I(2)*eta2;
cl2(:) = lhs2 - rhs2;

end
