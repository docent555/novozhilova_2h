function s = dfdt(t, f, u, n, Dtr, Ne, zax, Q, Th, A, Dr, R, I) %#codegen

s = zeros(6,1);

fp = complex(zeros(2,1));

fp(1) = f(1)*exp(1i*f(2));
fp(2) = f(3)*exp(1i*f(4));

Th0 = 2.0D0*pi*(0:Ne-1)/Ne;
p0 = exp(1i*Th0).';

opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

p1 = oscill_ode(fp(1), Ne, zax, Dtr(1), p0, u, n, opts);
p2 = oscill_ode(fp(2), Ne, zax, Dtr(2), p0, u, n, opts);

x1 = xi(p1, n, u, zax);
x2 = xi(p2, n, u, zax);

x1r = real(x1);
x1i = imag(x1);
x2r = real(x2);
x2i = imag(x2);

f1 = f(1);
phi1 = f(2);
f2 = f(3);
phi2 = f(4);
f3 = f(5);
phi3 = f(6);

Q31 = Q(3)/Q(1);
I1 = I(1);
R1 = R(1);
Th1 = Th(1);
Dr1 = Dr(1);
cos1 = cos(phi1);
sin1 = sin(phi1);

Q32 = Q(3)/Q(2);
I2 = I(2);
R2 = R(2);
Th2 = Th(2);
Dr2 = Dr(2);
cos2 = cos(phi2);
sin2 = sin(phi2);

Q3 = Q(3);
A1 = A(1);
A2 = A(2);

s(1) = (-n*f1 + I1*(-x1i*cos1 + x1r*sin1) + 2*R1*n*f3*cos(phi3 - phi1 - Th1))*Q31;
s(2) = -2*Dr1*Q3 + (I1/f1*(x1r*cos1 + x1i*sin1) + 2*R1*n*(f3/f1)*sin(phi3 - phi1 - Th1))*Q31;

s(3) = (-n*f2 + I2*(-x2i*cos2 + x2r*sin2) + 2*R2*n*f3*cos(phi3 - phi2 - Th2))*Q32;
s(4) = -2*Dr2*Q3 + (I2/f2*(x2r*cos2 + x2i*sin2) + 2*R2*n*(f3/f2)*sin(phi3 - phi2 - Th2))*Q32;

s(5) = -f3 + A1*f1*cos(phi1 - phi3) + A2*f2*cos(phi2 - phi3);
s(6) = A1*f1/f3*sin(phi1 - phi3) + A2*f2/f3*sin(phi2 - phi3);
end
