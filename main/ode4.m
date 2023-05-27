function data_out = ode4(odefun,tspan,y0,data)
%ODE4  Solve differential equations with a non-adaptive method of order 4.
%   Y = ODE4(ODEFUN,TSPAN,Y0) with TSPAN = [T1, T2, T3, ... TN] integrates
%   the system of differential equations y' = f(t,y) by stepping from T0 to
%   T1 to TN. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
%   The vector Y0 is the initial conditions at T0. Each row in the solution
%   array Y corresponds to a time specified in TSPAN.
%
%   Y = ODE4(ODEFUN,TSPAN,Y0,P1,P2...) passes the additional parameters
%   P1,P2... to the derivative function as ODEFUN(T,Y,P1,P2...).
%
%   This is a non-adaptive solver. The step sequence is determined by TSPAN
%   but the derivative function ODEFUN is evaluated multiple times per step.
%   The solver implements the classical Runge-Kutta method of order 4.
%
%   Example
%         tspan = 0:0.1:20;
%         y = ode4(@vdp1,tspan,[2 0]);
%         plot(tspan,y(:,1));
%     solves the system y' = vdp1(t,y) with a constant step size of 0.1,
%     and plots the first component of the solution.
silentMode=0; useWaitbar=1;
data_out=data;  %initialize output data
out_tspan=data.out_tspan;

%%
if ~isnumeric(tspan)
    error('TSPAN should be a vector of integration steps.');
end

if ~isnumeric(y0)
    error('Y0 should be a vector of initial conditions.');
end

h = diff(tspan);
if any(sign(h(1))*h <= 0)
    error('Entries of TSPAN are not in order.')
end

try
    f0 = feval(odefun,tspan(1),y0,data);
catch
    msg = ['Unable to evaluate the ODEFUN at t0,y0. ',lasterr];
    error(msg);
end

y0 = y0(:);   % Make a column vector.
if ~isequal(size(y0),size(f0))
    error('Inconsistent sizes of Y0 and f(t0,y0).');
end

neq = length(y0);
N = length(tspan);
Y = zeros(neq,N);
F = zeros(neq,4);

if silentMode==0 && useWaitbar==1
    wb= waitbar(0,'Please wait...');
end

%initialize time history data
data_out.Ya_t=zeros(neq,numel(out_tspan));
Y(:,1) = y0;
for i = 2:N
    ti = tspan(i-1);
    hi = h(i-1);
    yi = Y(:,i-1);
    F(:,1) = feval(odefun,ti,yi,data);
      if  sum(out_tspan==ti)
        if silentMode==0
            if useWaitbar==1
                waitbar(i/N,wb,['Current time: ',sprintf('%.5f',ti),'s']); % display a waitbar with the instant simulation clock
            else
                disp(['Current time: ',sprintf('%.5f',ti)]); % display the simulation clock
            end
            % disp(ti);
        end
        data_out.Ya_t(:,out_tspan==ti)=F(:,1);
       
    end
    F(:,2) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,1),data);
    F(:,3) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,2),data);
    F(:,4) = feval(odefun,tspan(i),yi+hi*F(:,3),data);
    Y(:,i) = yi + (hi/6)*(F(:,1) + 2*F(:,2) + 2*F(:,3) + F(:,4));
end
close(wb);
%% output data
% data_out.Ya_t=Y;
% data_out.na_t=Y(1:neq/2,:);
% data_out.na_d_t=Y(neq/2+1:end,:);

