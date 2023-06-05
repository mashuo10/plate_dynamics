%%% 加劲板动力时程分析_获取外力%%%%%%
%%%% 运行此程序，获取外力输入'PA_t_num','W_base','W_base_num','num_load'数据
clc;
clear all;
close all;



%% load K M n
load('K_M_matrix.mat');
% load('external_force_time.mat');
a=2.4;%x方向长度
b=0.4;%y方向长度

%荷载参数%%%
Cx=50;%荷载速度
tr=0.003;%荷载上升时间
td=0.0015;%荷载衰减时间
tc=tr+td;%荷载作用时间
lr=Cx*tr;%荷载上升宽度
ld=Cx*td;%荷载宽度
lc=Cx*tc;%荷载作用宽度
Pm=750000;%峰值荷载大小
Pl=75000;%低压荷载大小

%%
% time step
dt=1e-5;        %计算时间步长
out_dt=1e-5;        %输出时间步长
tf=0.1;
tspan=dt:dt:tf;
out_tspan=interp1(tspan,tspan,out_dt:out_dt:tf, 'nearest','extrap');  % output data time span


%% shape function
syms x y
ksi=2*x/a; yita=2*y/b;      % ksi and yita
% if ~(exist('W_base')  ) &&~(exist('W_base_num') )

%sin方函数构造
p1=[];
for i=1:n;
    p9=(-1)^i*(((ksi+1)/2)^3-((ksi+1)/2)^2)+((ksi+1)/2)^3-2*((ksi+1)/2)^2+(ksi+1)/2-1/i/pi*sin(i*pi*(ksi+1)/2);
    p1=[p1,p9];
end
tp1 = p1';

py1=[];
for i=1:n;
    p9=(-1)^i*(((yita+1)/2)^3-((yita+1)/2)^2)+((yita+1)/2)^3-2*((yita+1)/2)^2+(yita+1)/2-1/i/pi*sin(i*pi*(yita+1)/2);
    py1=[py1,p9];
end
tpy1 = py1.';

Int_y=integral(matlabFunction(tpy1),-b/2,b/2,'ArrayValued',true);   % integration of tpy1



% may delete this part
% W_base=cell(n,n);           % base matrix of W_ij
% W_base_num=cell(n,n);
% for i=1:n
%     for j=1:n
%         W_base{i,j}=((-1)^i*(((ksi+1)/2)^3-((ksi+1)/2)^2)+((ksi+1)/2)^3-2*((ksi+1)/2)^2+(ksi+1)/2-1/i/pi*sin(i*pi*(ksi+1)/2))*((-1)^j*(((yita+1)/2)^3-((yita+1)/2)^2)+((yita+1)/2)^3-2*((yita+1)/2)^2+(yita+1)/2-1/j/pi*sin(j*pi*(yita+1)/2));
%         W_base_num{i,j}=eval(subs(subs(W_base{i,j},x,linspace(-a/2,a/2,31)),y,linspace(-b/2,b/2,31)'));
%     end
% end

%% plot load P(x,y,t)
%{
if 0
    P_xyt=cell(numel(tspan),1);         %P(x,y,t)
    num=10;
plot_time=[0:5e-3:0.01,0.02:1e-2:0.12];
for k=1:num
    t=plot_time(k);                         % corresponding time
    P_xyt{k}=piecewise(Cx*t-lr-a/2<=x<Cx*t-a/2,Pm*(Cx*t-x-a/2)/lr,x<=Cx*t-lr-a/2,Pl+(Pm-Pl)*exp((x+a/2-Cx*t+lr)/ld),x>=Cx*t-a/2,0);  % define piecewise function分段荷载函数P(x,y,t)
    [X,Y] = meshgrid(linspace(-a/2,a/2,120),linspace(-b/2,b/2,20));
    Z=eval(subs(subs(P_xyt{k},x,linspace(-a/2,a/2,120)),y,linspace(-b/2,b/2,20)'));
    figure
    surf(X,Y,Z)
    xlabel('X(m)');
    ylabel('Y');
    zlabel('Pressure(Pa)')
    zlim([0,1.5*Pm]);
    title(['time=',num2str(t)]);
end
end
%}
%% calculat load coefficient PA_ij_t
dt_load=1e-4;                   % 每隔1e-4s算一次广义荷载
num_load=tf/dt_load;               % 
PA_t_num=zeros(n,n,num_load);    % PA_ij with different time in reduced number
for k=1:num_load
%     t=tspan(k*10);                         % corresponding time
    t=k*dt_load;
    P_xyt_line=(Pm*(Cx*t-x-a/2)/lr);                     %linear part
    P_xyt_exp=(Pl+(Pm-Pl)*exp((x+a/2-Cx*t+lr)/ld));     %exponential part
    if mod(k,1e2)==0
    fprintf('k =%d\n',k);
    end
    if t<=tr
        % aaa=eval(int(int(P_xyt{k}*W_base{i,j},x,-a/2,a/2),y,-b/2,b/2))
        % aaa=vpa(int(int(P_xyt{k}*W_base{i,j},x,-a/2,a/2),y,-b/2,b/2));      %符号积分
        % % 数值积分
        % fun = matlabFunction(P_xyt{k}*W_base{i,j});
        % aaa=integral2(fun,-a/2,a/2,-b/2,b/2,"AbsTol",1e-3)

        % double integral
%                         fun = matlabFunction(P_xyt_line*W_base{i,j});               %linear part
%                         aaa=integral2(fun,-a/2,Cx*t-a/2,-b/2,b/2,'AbsTol',1e-3);

        % integration twice based on linear independant of X Y
        % deformation
        Int_x=integral(matlabFunction(P_xyt_line*tp1),-a/2,Cx*t-a/2,'ArrayValued',true);
    elseif t>tr&& t<=a/Cx
        Int_x2=integral(matlabFunction(P_xyt_exp*tp1),-a/2,Cx*t-lr-a/2,'ArrayValued',true);%exponential part
        Int_x1=integral(matlabFunction(P_xyt_line*tp1),Cx*t-lr-a/2,Cx*t-a/2,'ArrayValued',true);%linear part
        Int_x=Int_x1+Int_x2;
    elseif t>a/Cx&& t<=a/Cx+tr
        Int_x2=integral(matlabFunction(P_xyt_exp*tp1),-a/2,Cx*t-lr-a/2,'ArrayValued',true); %exponential part
        Int_x1=integral(matlabFunction(P_xyt_line*tp1),Cx*t-lr-a/2,a/2,'ArrayValued',true);%linear part
        Int_x=Int_x1+Int_x2;
    elseif t>a/Cx+tr
        Int_x2=integral(matlabFunction(P_xyt_exp*tp1),-a/2,a/2,'ArrayValued',true);%exponential part
        Int_x=Int_x2;
    end
    PA_t_num(:,:,k)=Int_x*Int_y';
end

PA_t=zeros(n,n,numel(tspan));    % PA_ij with different time in real time
for i=1:numel(tspan)
%     PA_t(:,:,i)=PA_t_num(:,:,ceil((i/numel(tspan))*num_load));
    PA_t(:,:,i)=PA_t_num(:,:,ceil(i*dt/dt_load));
end

save (['external_force_time','.mat'],'PA_t_num','W_base','W_base_num','num_load');

%% plot load P(x,y,t)   (this part may be wrong)
plot_num=round(linspace(1,numel(tspan),8));

% calculate displacement
Z_disp=cell(numel(plot_num),1);                    % displacement of whole plate
for k=1:numel(plot_num)
    Z_disp{k}=zeros(31,31);
    
    for i=1:n
        for j=1:n
            Z_disp{k}=Z_disp{k}+W_base_num{i,j}*PA_t(i,j,plot_num(k));%displacement of whole plate
            
        end
    end
end

[X,Y] = meshgrid(linspace(-a/2,a/2,31),linspace(-b/2,b/2,31));

for k=1:numel(plot_num)
% [~,k]=find(tspan==plot_time(i));
% k=plot_num(i);
figure 
    surf(X,Y,Z_disp{k})
    xlabel('X(m)');
    ylabel('Y(m)');
    zlabel('Z坐标(m)')
%     zlim([0,1.5*Pm]);
    title(['time=',num2str(tspan(plot_num(k)))]);
end
%}


%% static analysis

K_num_mod=K_num*mod_K;

M_num_mod=M_num*mod_M;

%%  calculate the displacement in evenly distributed load
Int_x=integral(matlabFunction(Pm*tp1),-a/2,a/2,'ArrayValued',true);

PA_junbu=Int_x*Int_y';      % load factor in evenly distributed load
Z_disp_vect=K_num_mod\reshape(PA_junbu',n*n,1);

Z_disp_junbu=zeros(31,31);
    for i=1:n
        for j=1:n
           Z_disp_junbu=Z_disp_junbu+W_base_num{i,j}*Z_disp_vect((i-1)*n+j);%displacement of whole plate
            
        end
    end

figure
 [X,Y] = meshgrid(linspace(-a/2,a/2,31),linspace(-b/2,b/2,31));
    surf(X,Y,Z_disp_junbu)
    xlabel('X(m)');
    ylabel('Y(m)');
    zlabel('Z坐标(m)');
    title('均布荷载位移');
%     zlim([0,1.5*Pm]);

max(Z_disp_junbu(:))



 
%% calculate dynamics response

% initial value
a0=zeros(n*n,1);    % initial displacement  a0
a0_d=zeros(n*n,1);  % initial velocity a0_d
Y0a=[a0;a0_d];      %initial value of a0, a0_d

%damping matrix
D=0*eye(n*n);

%% data_in

% % redistribute PA_t from PA_t_num
% PA_t=zeros(n,n,numel(tspan));    % PA_ij with different time in real time
% for i=1:numel(tspan)
%     PA_t(:,:,i)=PA_t_num(:,:,ceil(i/(numel(tspan)/num_load)));
% end

% external force
w_t=zeros(n*n,numel(tspan));
for i=1:numel(tspan)
    w_t(:,i)=reshape(PA_t(:,:,i)',n*n,1);        % be careful of the row and column order!!!!!!!!!
%     w_t(:,i)=i/numel(tspan)*reshape(PA_junbu',n*n,1);              % use evenly distributed load    
end

global data

data.tf=tf;data.dt=dt;data.tspan=tspan;data.out_tspan=out_tspan;
data.M=M_num_mod;data.K=K_num_mod; data.D=D;
data.w_t=w_t;

%% solve dynamics equations
tic;

%  data_out = ode4(@Y_dot,tspan,Y0a,data);
%  Ya_t=data_out.Ya_t;

 [t,y]=ode23(@Y_dot_ty,tspan,Y0a);
Ya_t=y';

toc;
nf=n*n;
na_t=Ya_t(1:nf,:);               %free node cooridnate
na_d_t=Ya_t(nf+1:end,:);         %free node velocity
%% save data
save(['dynamics_result_ode23_0.1s_',num2str(1e3*a),'_',num2str(1e3*b),'.mat'],'na_t','na_d_t','out_tspan');
%% calculate displacement
Z_disp_x_0=zeros(31,numel(out_tspan));           % displacement in x=0
Z_disp_y_0=zeros(31,numel(out_tspan));           % displacement in y=0
Z_disp=cell(numel(out_tspan),1);                    % displacement of whole plate
for k=1:numel(out_tspan)
    Z_disp{k}=zeros(31,31);
    for i=1:n
        for j=1:n
            Z_disp{k}=Z_disp{k}+W_base_num{i,j}*na_t((i-1)*n+j,k);%displacement of whole plate
            
        end
    end
Z_disp_x_0(:,k)=Z_disp{k}(:,16) ; % displacement in x=0
Z_disp_y_0(:,k)=Z_disp{k}(16,:)' ; % displacement in y=0
end

%% plot displacement 
% center node displacement
figure
plot(out_tspan,Z_disp_x_0(16,:),'LineWidth',1.5);
grid on
xlabel('时间(s)');
ylabel('位移(m)');
title('中心点位移时程曲线');

fprintf(['maximum displacement is: ',num2str(1e3*max(Z_disp_x_0(16,:))),'mm'])

return



% X=0 displcement

figure
[X,Y] = meshgrid(out_tspan,linspace(-b/2,b/2,31));
s=surf(X,Y,Z_disp_x_0);
s.EdgeColor = 'none';
xlabel('时间(s)');
ylabel('Y坐标(m)');
    zlabel('Z坐标(m)')
title('X=0位移时程曲线');

% Y=0 displcement

figure
[X,Y] = meshgrid(out_tspan,linspace(-a/2,a/2,31));
s=surf(X,Y,Z_disp_y_0);
s.EdgeColor = 'none';
xlabel('时间(s)');
ylabel('X坐标(m)');
    zlabel('Z坐标(m)')
title('Y=0位移时程曲线');

%% plot plate configuration
plot_num=round(linspace(1,numel(out_tspan),8));

for i=1:numel(plot_num)
% [~,k]=find(tspan==plot_time(i));
k=plot_num(i);
figure
 [X,Y] = meshgrid(linspace(-a/2,a/2,31),linspace(-b/2,b/2,31));
    surf(X,Y,Z_disp{k})
    xlabel('X(m)');
    ylabel('Y(m)');
    zlabel('Z坐标(m)')
%     zlim([0,1.5*Pm]);
    title(['time=',num2str(out_tspan(k))]);
end
%}




























return
%%
plot_time=[2e-3:2e-3:0.01,0.02:1e-2:0.12];
Z_disp=cell(numel(out_tspan),1);                    % displacement of whole plate
for i=1:numel(plot_time)
[~,k]=find(tspan==plot_time(i));
Z_disp_vect=K_num_mod\reshape(PA_t(:,:,k)',n*n,1);
    Z_disp{k}=zeros(31,31);
    for i=1:n
        for j=1:n
           Z_disp{k}=Z_disp{k}+W_base_num{i,j}*Z_disp_vect((i-1)*n+j,k);%displacement of whole plate
            
        end
    end

figure
 [X,Y] = meshgrid(linspace(-a/2,a/2,31),linspace(-b/2,b/2,31));
    surf(X,Y,Z_disp{k})
    xlabel('X(m)');
    ylabel('Y(m)');
    zlabel('Z坐标(m)')
    zlim([0,1.5*Pm]);
    title(['time=',num2str(plot_time(i))]);
end


