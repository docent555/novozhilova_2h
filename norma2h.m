clear all

indata=read_namelist('input_norma_gmp.in','PARAM');
pitch   =  indata.pitch;
nharm   =  indata.nharm;
UKV     =  indata.ukv;
Fop     =  indata.f_op;
m       =  indata.m;
q       =  indata.q;
Rb      =  indata.rb;
Rr      =  indata.rr;
estx    =  indata.estx;

gamma = 1 + UKV/511;
Beta = sqrt(1 - 1/(gamma*gamma));
Beta2 = 1 - 1/(gamma*gamma);
Beta_z = Beta/sqrt(pitch*pitch + 1);
Beta_z2 = Beta2/(pitch*pitch + 1.0d0);
Beta_perp2 = Beta2 - Beta_z2;
Beta_perp = sqrt(Beta_perp2);
Wop = 2*pi*Fop*1e9;
c = 299792458*10^2;

struct = load('ColdTabl_28GHz_bezzgolovka.dat');

u = struct(:,3).*exp(1i*struct(:,4));
dz = 0.0280211;
dz = Beta_perp2/2/Beta_z*Wop*dz/nharm/c;
norm = sum(abs(u(:)).*abs(u(:)))*dz

Gmp = gmp(m, q, Rb, Rr, estx)