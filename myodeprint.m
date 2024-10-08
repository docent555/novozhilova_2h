function status = myodeprint(t,y,flag,varargin) %#codegen
%ODEPRINT  Command window printing ODE output function.
%   When the function odeprint is passed to an ODE solver as the 'OutputFcn'
%   property, i.e. options = odeset('OutputFcn',@odeprint), the solver calls
%   ODEPRINT(T,Y,'') after every timestep. The ODEPRINT function prints all
%   components of the solution it is passed as it is computed. To print only
%   particular components, specify their indices in the 'OutputSel' property
%   passed to the ODE solver.
%
%   At the start of integration, a solver calls ODEPRINT(TSPAN,Y0,'init') to
%   initialize the output function.  After each integration step to new time
%   point T with solution vector Y the solver calls STATUS = ODEPRINT(T,Y,'').
%   If the solver's 'Refine' property is greater than one (see ODESET), then
%   T is a column vector containing all new output times and Y is an array
%   comprised of corresponding column vectors.  ODEPRINT always returns
%   STATUS = 0.  When the integration is complete, the solver calls
%   ODEPRINT([],[],'done').
%
%   See also ODEPLOT, ODEPHAS2, ODEPHAS3, ODE45, ODE15S, ODESET.

%   Mark W. Reichelt and Lawrence F. Shampine, 3-24-94
%   Copyright 1984-2020 The MathWorks, Inc.

persistent fidF fidP fidE fidCl Tend Zex Ne Nz Nt Q ...
    I Th A Dr Dtr R dt dz nharm zax u fp th0 p0 opts eta1 eta2

if isempty(fidF)
    fidF = double(1);
    fidP = double(1);    
    fidE = double(1);
    fidCl = double(1);   

    opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

    fp = complex(zeros(2,1));
    Ne = zeros(1,1);
    Zex = zeros(1,1);
    Dtr = zeros(1,1);    
    p0 = zeros(1,1);
    u = griddedInterpolant((1:663)',complex((1:663))');
    nharm = zeros(1,1);
    R = zeros(2,1);
    Th = zeros(2,1);
    I = zeros(2,1);
    eta1 = zeros(1,1);
    eta2 = zeros(1,1);
end

if nargin < 3 || isempty(flag) % odeprint(t,y) [v5 syntax] or odeprint(t,y,'')

    % clc
    % if ~isempty(t)
    %     t(end)   %#ok<NOPRT>
    % end
    % if ~isempty(y)
    %     y(1,end)  %#ok<NOPRT>
    % end
    if ~isempty(t) && ~isempty(y)

        fprintf(['\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
            '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
            't = %8.5f |F1| = %+8.5f |F2| = %+8.5f |F3| = %+8.5f'], t(end), y(1,end), y(3,end), y(5,end));

        fprintf(fidF,'%17.8e\t%17.8e\t%17.8e\t%17.8e\n',  t(end), y(1,end), y(3,end), y(5,end));
        fprintf(fidP,'%17.8e\t%17.8e\t%17.8e\t%17.8e\n',  t(end), y(2,end), y(4,end), y(6,end));

        n = length(t);

        % Вычисление КПД
        for i=1:n
            fp(1) = y(1,i)*exp(1i*y(2,i));
            fp(2) = y(3,i)*exp(1i*y(4,i));

            p1 = oscill_ode(fp(1), Ne, [0 Zex], Dtr(1), p0, u, nharm, opts);
            p2 = oscill_ode(fp(2), Ne, [0 Zex], Dtr(2), p0, u, nharm, opts);

            eta1 = 1 - mean(abs(p1(end,:)).^2);
            eta2 = 1 - mean(abs(p2(end,:)).^2);

            fprintf(fidE,'%17.8e\t%17.8e\t%17.8e\n', t(i), eta1, eta2);
        end
        %-------------------------------

        % Баланс
        for i=1:n
            lhs1 = 2*nharm*y(1,i).^2 - 4*nharm*R(1)*y(5,i).*y(1,i).*cos(Th(1) - y(6,i) + y(2,i));
            rhs1 = I(1)*eta1;
            cl1 = lhs1 - rhs1;

            lhs2 = 2*nharm*y(3,i).^2 - 4*nharm*R(2)*y(5,i).*y(3,i).*cos(Th(2) - y(6,i) + y(4,i));
            rhs2 = I(2)*eta2;
            cl2 = lhs2 - rhs2;

            fprintf(fidCl,'%17.8e\t%17.8e\t%17.8e\n', t(i), cl1, cl2);
        end
    end

else
    if isstring(flag) && isscalar(flag)
        flag = char(flag);
    end
    switch(flag)
        case 'init'               % odeprint(tspan,y0,'init')

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

            % clc
            % t = t(1) %#ok<NASGU,NOPRT>
            % y  %#ok<NOPRT>

            fidF = fopen('F.dat','w');
            fidP = fopen('P.dat','w');

            % Для КПД
            fidE = fopen('E.dat','w');            
            th0 = 2.0D0*pi*(0:Ne-1)/Ne;
            p0 = exp(1i*th0).';
            % opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

            fp(1) = y(1,1)*exp(1i*y(2,1));
            fp(2) = y(3,1)*exp(1i*y(4,1));

            p1 = oscill_ode(fp(1), Ne, [0 Zex], Dtr(1), p0, u, nharm, opts);
            p2 = oscill_ode(fp(2), Ne, [0 Zex], Dtr(2), p0, u, nharm, opts);

            eta1 = 1 - mean(abs(p1(end,:)).^2);
            eta2 = 1 - mean(abs(p2(end,:)).^2);

            fprintf(fidE,'%17.8e\t%17.8e\t%17.8e\n', t(1), eta1, eta2);
            %-----------------------------------

            % Баланс
            fidCl = fopen('cl.dat','w');

            lhs1 = 2*nharm*y(1,1).^2 - 4*nharm*R(1)*y(5,1).*y(1,1).*cos(Th(1) - y(6,1) + y(2,1));
            rhs1 = I(1)*eta1;
            cl1 = lhs1 - rhs1;

            lhs2 = 2*nharm*y(3,1).^2 - 4*nharm*R(2)*y(5,1).*y(3,1).*cos(Th(2) - y(6,1) + y(4,1));
            rhs2 = I(2)*eta2;
            cl2 = lhs2 - rhs2;

            fprintf(fidCl,'%17.8e\t%17.8e\t%17.8e\n', t(1), cl1, cl2);
            %-----------------------------------
        case 'done'               % odeprint([],[],'done')

            % fprintf('\n\n');

            fclose(fidF);
            fclose(fidP);
            fclose(fidE);
            fclose(fidCl);

    end
end

status = 0;
end
