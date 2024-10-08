function solve_odeF(Tend, Zex, Ne, Nz, Q, ...
    I, Th, A, Dr, Dtr, R, F0, P0, dt, dz, nharm) %#codegen

Nt = Tend/dt + 1;
if Nz > 0
    dz = Zex/(Nz - 1);
    % disp('nz > 0')
else
    Nz = fix(Zex/dz) + 1;
    % disp('nz == 0')
end

tax = 0:dt:Tend;
zax = 0:dz:Zex;

% eta1 = zeros(Nt,1);
% eta2 = zeros(Nt,1);
% cl1 = zeros(Nt,1);
% cl2 = zeros(Nt,1);

u = calc_u(Zex);

% параметры для расчета якобиана
fu = u(0:dz:Zex);
fu = [(0:dz:Zex)' real(fu)' imag(fu)'];

fid = fopen('u.bin','w');
fwrite(fid,fu,'double');
fclose(fid);

par = [Tend, Zex, Ne, Nz, Nt, Q, ...
    I, Th, A, Dr, Dtr, R, dt, dz, nharm];

fid = fopen('par.bin','w');
fwrite(fid,par,'double');
fclose(fid);
%----------------------------------

f0 = [F0(1); P0(1); F0(2); P0(2); F0(3); P0(3)];

opts = odeset('RelTol',1e-5,'AbsTol',1e-8,'OutputFcn',@myodeprint,'Jacobian',@jacobian);

fprintf('t = %8.5f |F1| = %+8.5f |F2| = %+8.5f |F3| = %+8.5f', 0, f0(1), f0(3), f0(5));
[~, f] = ode15s(@(t,f) dfdt(t, f, u, nharm, Dtr, Ne, zax, Q, Th, A, Dr, R, I), tax, f0, opts);
fprintf('\n');

% [eta1, eta2] = eff(eta1, eta2, f, Dtr, Ne, Nt, Zex, u, nharm);
% [cl1, cl2] = clow(cl1, cl2, eta1, eta2, f, Th, R, nharm, I);

% fileID = fopen('F.dat','w');
% for i=1:Nt
%     fprintf(fileID,'%17.8e\t%17.8e\t%17.8e\t%17.8e\n', tax(i), f(i,1), f(i,3), f(i,5)); 
% end
% fclose(fileID);

% fileID = fopen('P.dat','w');
% for i=1:Nt
%     fprintf(fileID,'%17.8e\t%17.8e\t%17.8e\t%17.8e\n', tax(i), f(i,2), f(i,4), f(i,6)); 
% end
% fclose(fileID);

% fileID = fopen('E1.dat','w');
% for i=1:Nt
%     fprintf(fileID,'%17.8e\t%17.8e\t%17.8e\n', tax(i), eta1(i), eta2(i)); 
% end
% fclose(fileID);

% fileID = fopen('cl1.dat','w');
% for i=1:Nt
%     fprintf(fileID,'%17.8e\t%17.8e\t%17.8e\n', tax(i), cl1(i), cl2(i)); 
% end
% fclose(fileID);
% end
