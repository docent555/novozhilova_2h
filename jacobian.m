function dfdy = jacobian(t,y) %#codegen

persistent Tend Zex Ne Nz Nt Q ...
    I Th A Dr Dtr R dt dz nharm zax u Q31 I1 R1 Th1 A1 Q32 I2 R2 Th2 A2

if isempty(nharm)
    Q = zeros(3,1);
    I = zeros(2,1);
    Th = zeros(2,1);
    A = zeros(2,1);
    Dr = zeros(2,1);
    Dtr = zeros(2,1);
    R = zeros(2,1);

    fid = fopen('par.bin');
    par = fread(fid, [23 1],'double');
    fclose(fid);
    
    Tend = par(1);
    Zex = par(2);
    Ne = par(3);
    Nz = par(4);
    Nt = par(5);
    Q(1) = par(6);
    Q(2) = par(7);
    Q(3) = par(8);
    I(1) = par(9);
    I(2) = par(10);
    Th(1) = par(11);
    Th(2) = par(12);   
    A(1) = par(13);
    A(2) = par(14);
    Dr(1) = par(15);
    Dr(2) = par(16);
    Dtr(1) = par(17);
    Dtr(2) = par(18);
    R(1) = par(19);
    R(2) = par(20);    
    dt = par(21);
    dz = par(22);
    nharm = par(23);

    zax = 0:dz:Zex;  
    u = calc_u(Zex);

    Q31 = Q(3)/Q(1);
    I1 = I(1);
    R1 = R(1);
    Th1 = Th(1);

    Q32 = Q(3)/Q(2);
    I2 = I(2);
    R2 = R(2);
    Th2 = Th(2);

    % Q3 = Q(3);
    A1 = A(1);
    A2 = A(2);
end

dfdy = zeros(6);

fp = complex(zeros(2,1));

fp(1) = y(1)*exp(1i*y(2));
fp(2) = y(3)*exp(1i*y(4));

Th0 = 2.0D0*pi*(0:Ne-1)/Ne;
p0 = exp(1i*Th0).';

opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

p1 = oscill_ode(fp(1), Ne, zax, Dtr(1), p0, u, nharm, opts);
p2 = oscill_ode(fp(2), Ne, zax, Dtr(2), p0, u, nharm, opts);

x1 = xi(p1, nharm, u, zax);
x2 = xi(p2, nharm, u, zax);

x1r = real(x1);
x1i = imag(x1);
x2r = real(x2);
x2i = imag(x2);

f1 = y(1);
phi1 = y(2);
f2 = y(3);
phi2 = y(4);
f3 = y(5);
phi3 = y(6);

% Q31 = Q(3)/Q(1);
% I1 = I(1);
% R1 = R(1);
% Th1 = Th(1);
% 
% Q32 = Q(3)/Q(2);
% I2 = I(2);
% R2 = R(2);
% Th2 = Th(2);
% 
% % Q3 = Q(3);
% A1 = A(1);
% A2 = A(2);

dfdy(1,1) = (-1).*nharm.*Q31;
dfdy(1,2) = I1.*Q31.*x1r.*cos(phi1)+I1.*Q31.*x1i.*sin(phi1)+(-2).*f3.*nharm.* ...
  Q31.*R1.*sin(phi1+(-1).*phi3+Th1);
dfdy(1,3) = 0;
dfdy(1,4) = 0;
dfdy(1,5) = 2.*nharm.*Q31.*R1.*cos(phi1+(-1).*phi3+Th1);
dfdy(1,6) = 2.*f3.*nharm.*Q31.*R1.*sin(phi1+(-1).*phi3+Th1);

dfdy(2,1) = (-1).*f1.^(-2).*I1.*Q31.*x1r.*cos(phi1)+(-1).*f1.^(-2).*I1.*Q31.* ...
  x1i.*sin(phi1)+2.*f1.^(-2).*f3.*nharm.*Q31.*R1.*sin(phi1+(-1).* ...
  phi3+Th1);
dfdy(2,2) = f1.^(-1).*I1.*Q31.*x1i.*cos(phi1)+(-2).*f1.^(-1).*f3.*nharm.*Q31.* ...
  R1.*cos(phi1+(-1).*phi3+Th1)+(-1).*f1.^(-1).*I1.*Q31.*x1r.*sin( ...
  phi1);
dfdy(2,3) = 0;
dfdy(2,4) = 0;
dfdy(2,5) = (-2).*f1.^(-1).*nharm.*Q31.*R1.*sin(phi1+(-1).*phi3+Th1);
dfdy(2,6) = 2.*f1.^(-1).*f3.*nharm.*Q31.*R1.*cos(phi1+(-1).*phi3+Th1);

dfdy(3,1) = 0;
dfdy(3,2) = 0;
dfdy(3,3) = (-1).*nharm.*Q32;
dfdy(3,4) = I2.*Q32.*x2r.*cos(phi2)+I2.*Q32.*x2i.*sin(phi2)+(-2).*f3.*nharm.* ...
  Q32.*R2.*sin(phi2+(-1).*phi3+Th2);
dfdy(3,5) = 2.*nharm.*Q32.*R2.*cos(phi2+(-1).*phi3+Th2);
dfdy(3,6) = 2.*f3.*nharm.*Q32.*R2.*sin(phi2+(-1).*phi3+Th2);

dfdy(4,1) = 0;
dfdy(4,2) = 0;
dfdy(4,3) = (-1).*f2.^(-2).*I2.*Q32.*x2r.*cos(phi2)+(-1).*f2.^(-2).*I2.*Q32.* ...
  x2i.*sin(phi2)+2.*f2.^(-2).*f3.*nharm.*Q32.*R2.*sin(phi2+(-1).* ...
  phi3+Th2);
dfdy(4,4) = f2.^(-1).*I2.*Q32.*x2i.*cos(phi2)+(-2).*f2.^(-1).*f3.*nharm.*Q32.* ...
  R2.*cos(phi2+(-1).*phi3+Th2)+(-1).*f2.^(-1).*I2.*Q32.*x2r.*sin( ...
  phi2);
dfdy(4,5) = (-2).*f2.^(-1).*nharm.*Q32.*R2.*sin(phi2+(-1).*phi3+Th2);
dfdy(4,6) = 2.*f2.^(-1).*f3.*nharm.*Q32.*R2.*cos(phi2+(-1).*phi3+Th2);

dfdy(5,1) = A1.*cos(phi1+(-1).*phi3);
dfdy(5,2) = (-1).*A1.*f1.*sin(phi1+(-1).*phi3);
dfdy(5,3) = A2.*cos(phi2+(-1).*phi3);
dfdy(5,4) = (-1).*A2.*f2.*sin(phi2+(-1).*phi3);
dfdy(5,5) = (-1);
dfdy(5,6) = A1.*f1.*sin(phi1+(-1).*phi3)+A2.*f2.*sin(phi2+(-1).*phi3);

dfdy(6,1) = A1.*f3.^(-1).*sin(phi1+(-1).*phi3);
dfdy(6,2) = A1.*f1.*f3.^(-1).*cos(phi1+(-1).*phi3);
dfdy(6,3) = A2.*f3.^(-1).*sin(phi2+(-1).*phi3);
dfdy(6,4) = A2.*f2.*f3.^(-1).*cos(phi2+(-1).*phi3);
dfdy(6,5) = (-1).*A1.*f1.*f3.^(-2).*sin(phi1+(-1).*phi3)+(-1).*A2.*f2.*f3.^( ...
  -2).*sin(phi2+(-1).*phi3);
dfdy(6,6) = (-1).*A1.*f1.*f3.^(-1).*cos(phi1+(-1).*phi3)+(-1).*A2.*f2.*f3.^( ...
  -1).*cos(phi2+(-1).*phi3);

end

