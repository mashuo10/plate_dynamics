%%% 加劲板动力时程分析_获取外力%%%%%%
%%%% 运行此程序，获取外力输入'PA_t_num','W_base','W_base_num','num_load'数据
clc;
clear all;
close all;



%% load K M n
load('K_M_matrix.mat');
load('external_force_time.mat');
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
dt=1e-4;
out_dt=1e-3;
tf=0.04;
tspan=dt:dt:tf;
out_tspan=interp1(tspan,tspan,0:out_dt:tf, 'nearest','extrap');  % output data time span



%% shape function
syms x y
%if ~(exist('W_base')  ) &&~(exist('W_base_num') )

W_base=cell(n,n);           % base matrix of W_ij
W_base_num=cell(n,n);
ksi=2*x/a; yita=2*y/b;      % ksi and yita
for i=1:n
    for j=1:n
        W_base{i,j}=((-1)^i*(((ksi+1)/2)^3-((ksi+1)/2)^2)+((ksi+1)/2)^3-2*((ksi+1)/2)^2+(ksi+1)/2-1/i/pi*sin(i*pi*(ksi+1)/2))*((-1)^j*(((yita+1)/2)^3-((yita+1)/2)^2)+((yita+1)/2)^3-2*((yita+1)/2)^2+(yita+1)/2-1/j/pi*sin(j*pi*(yita+1)/2));
        W_base_num{i,j}=eval(subs(subs(W_base{i,j},x,linspace(-a/2,a/2,31)),y,linspace(-b/2,b/2,31)'));
    end
end



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
num_load=tf/1e-3;               % 每隔1e-3s算一次广义荷载
PA_t_num=zeros(n,n,num_load);    % PA_ij with different time in reduced number
for k=1:num_load
    t=tspan(k*10);                         % corresponding time
    P_xyt_line=Pm;                    %linear part
    for i=1:n
        fprintf('k =%d i= %d \n',k,i);
        for j=1:n
            if t<=tr
                % aaa=eval(int(int(P_xyt{k}*W_base{i,j},x,-a/2,a/2),y,-b/2,b/2))
                % aaa=vpa(int(int(P_xyt{k}*W_base{i,j},x,-a/2,a/2),y,-b/2,b/2));      %符号积分
                % % 数值积分
                % fun = matlabFunction(P_xyt{k}*W_base{i,j});
                % aaa=integral2(fun,-a/2,a/2,-b/2,b/2,"AbsTol",1e-3)

                fun = matlabFunction(P_xyt_line*W_base{i,j});               %linear part
                aaa=integral2(fun,-a/2,Cx*t-a/2,-b/2,b/2,'AbsTol',1e-3);
            elseif t>tr&& t<=a/Cx
                fun = matlabFunction(P_xyt_line*W_base{i,j});           %linear part
                aaa1=integral2(fun,Cx*t-lr-a/2,Cx*t-a/2,-b/2,b/2,'AbsTol',1e-3);
                aaa=aaa1;
            elseif t>a/Cx&& t<=a/Cx+tr
                fun = matlabFunction(P_xyt_line*W_base{i,j});           %linear part
                aaa1=integral2(fun,Cx*t-lr-a/2,a/2,-b/2,b/2,'AbsTol',1e-3);
                aaa=aaa1;
            end
            PA_t_num(i,j,k)=aaa;
            % integral2(matlabFunction(ttt),-a/2,a/2,-b/2,b/2);     %this doesn't work
        end
    end
end

PA_t=zeros(n,n,numel(tspan));    % PA_ij with different time in real time
for i=1:numel(tspan)
    PA_t(:,:,i)=PA_t_num(:,:,ceil(i/(numel(tspan)/num_load)));
end

save (['external_force_time','.mat'],'PA_t_num','W_base','W_base_num','num_load');



return
%% static analysis
%{
plot_time=[2e-3:2e-3:0.01,0.02:1e-2:0.12];
Z_disp=cell(numel(out_tspan),1);                    % displacement of whole plate
for i=1:numel(plot_time)
[~,k]=find(tspan==plot_time(i));
Z_disp_vect=K_num\reshape(PA_t(:,:,k)',n*n,1);
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
%     zlim([0,1.5*Pm]);
    title(['time=',num2str(plot_time(i))]);
end
%}

