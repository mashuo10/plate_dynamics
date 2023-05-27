function Yd=Y_dot_ty(t,Y)
% Input:
%   t: current time
%   Y0: current X, Xd values:   Y0=[X;Xd];
% Output:
%   Yd=[Xd,Xdd]
global data

w_t=data.w_t;         % external force
K=data.K;             %stiffness matrix
M=data.M;             % mass matrix
D=data.D;             %damping matrix
dt=data.dt;

nf=numel(Y)/2;
na=Y(1:nf,:);               %free node cooridnate
na_d=Y(nf+1:end,:);         %free node velocity
%% external force
ind=floor(t/dt);

% Get current external forces
if ischar(w_t)
    run(w_t) % if w_t is a string, run that script
elseif size(w_t,2)==1
    w = w_t; % w_t can be constant
else
    w = w_t(:,ind); % or W can be time-varying
end

%% calculate accerlation
na_dd=M\(w-K*na-D*na_d);      %dynamic equation
Yd=[na_d;na_dd];
