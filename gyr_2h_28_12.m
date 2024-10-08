clearvars

indata=read_namelist('input_fortran.in','PARAM');
Ne    = indata.ne;
Tend  = indata.tend;
% Zex   = indata.zex;
Q     = indata.q;
I     = indata.i;
th   = indata.th;
A     = indata.a;
Dr    = indata.dr;
R     = indata.r;
F0    = indata.f0;
PHI0  = indata.p0;
dt    = indata.dt;
dz    = indata.dz;
Dtrb  = indata.dtrb;
Dtrh  = indata.dtrh;
inher = indata.inher;
Nz    = indata.nz;

nharm = 2;
Zex = 9.6638;

tic
for dtr = Dtrb(1):Dtrh:Dtrb(2)
    Dtr = [dtr dtr];
    fprintf('Dtr1 = %f\tDtr2 = %f\n', dtr, dtr);
    solve_odeF_mex(Tend, Zex, Ne, Nz, Q, I, th, A, Dr, ...
        Dtr, R, F0, PHI0, dt, dz, nharm);
end
toc