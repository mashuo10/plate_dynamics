%%% 加劲板动力时程分析_获取质量矩阵，刚度矩阵%%%%%%
%%%% 运行此程序，获取刚度矩阵、质量矩阵、修正系数'K_num','M_num','n','mod_K','mod_M'数据
clc ;
clear all;
close all;

tic
%变量定义
syms z o
syms x y real
%阶数输入
n=15;
%板参数输入
a=1.5;%x方向长度
b=1.5;%y方向长度
h=0.015;%板厚
E=2.1e11;%弹性模量
v=0.3;%泊松比
P=7850;%材料密度
%x方向梁参数输入
Hx=0.14555;%腹板高度
Bx=0.0485;%翼缘宽度
t1x=0.0065;%腹板厚度
t2x=0.0145;%翼缘厚度
yk=[-0.25 0.25];%x方向梁位置
k=numel(yk);%x方向梁数量
% Hx=0;%腹板高度
% Bx=0;%翼缘宽度
% t1x=0;%腹板厚度
% t2x=0;%翼缘厚度
% yk=[];%x方向梁位置
% k=numel(yk);%x方向梁数量
%y方向梁参数输入
Hy=0;%腹板高度
By=0;%翼缘宽度
t1y=0;%腹板厚度
t2y=0;%翼缘厚度
xl=[];%y方向梁位置
l=numel(xl);%y方向梁数量

G=E/(2*(1+v));
D=E*h^3/(12*(1-v^2));  
r=a/b;

Ix=1/12*t1x*Hx^3+Hx*t1x*(Hx/2+h/2)^2+1/12*Bx*t2x^3+t2x*Bx*(Hx+t2x/2+h/2)^2;%center of top plates
% Ix=1/12*t1x*Hx^3+1/12*Bx*t2x^3+t2x*Bx*(Hx/2+t2x/2)^2;     %center of vertical plate


Ix1=1/12*(Hx*t1x^3+t2x*Bx^3);
% Jx=Ix+Ix1;  
Jx=1/3*(Hx*t1x^3+Bx*t2x^3);
Ax=t1x*Hx+t2x*Bx;
yk=yk*2/b;

Iy=1/12*t1y*Hy^3+Hy*t1y*(Hy/2+h/2)^2+1/12*By*t2y^3+t2y*By*(Hy+t2y/2+h/2)^2;
Iy1=1/12*(Hy*t1y^3+t2y*By^3);
% Jy=Iy+Iy1;  
Jy=1/3*(Hy*t1y^3+By*t2y^3);
Ay=t1y*Hy+t2y*By;
xl=xl*2/a;

%sin方函数构造
p1=[];
for i=1:n;
    p9=(-1)^i*(((x+1)/2)^3-((x+1)/2)^2)+((x+1)/2)^3-2*((x+1)/2)^2+(x+1)/2-1/i/pi*sin(i*pi*(x+1)/2);
    p1=[p1,p9];
end
tp1 = p1';

py1=[];
for i=1:n;
    p9=(-1)^i*(((y+1)/2)^3-((y+1)/2)^2)+((y+1)/2)^3-2*((y+1)/2)^2+(y+1)/2-1/i/pi*sin(i*pi*(y+1)/2);
    py1=[py1,p9];
end
tpy1 = py1.';

p00=[];
for i=1:n;
    p9=(-1)^i*(((z+1)/2)^3-((z+1)/2)^2)+((z+1)/2)^3-2*((z+1)/2)^2+(z+1)/2-1/i/pi*sin(i*pi*(z+1)/2);
    p00=[p00,p9];
end
tp00 = p00.';

%kesi
% Ft00t=int(tp1*p1,x,-1,1);
% Ft11t=int(diff(tp1,x)*diff(p1,x),x,-1,1);
% Ft02t=int(tp1*diff(p1,x,2),x,-1,1);
% Ft20t=int(diff(tp1,x,2)*p1,x,-1,1);
% Ft22t=int(diff(tp1,x,2)*diff(p1,x,2),x,-1,1);
% 
Ft00t=integral(matlabFunction(tp1*p1),-1,1,'ArrayValued',true);       %numerical integration is much faster
Ft11t=integral(matlabFunction(diff(tp1,x)*diff(p1,x)),-1,1,'ArrayValued',true);
Ft02t=integral(matlabFunction(tp1*diff(p1,x,2)),-1,1,'ArrayValued',true);
Ft20t=integral(matlabFunction(diff(tp1,x,2)*p1),-1,1,'ArrayValued',true);
Ft22t=integral(matlabFunction(diff(tp1,x,2)*diff(p1,x,2)),-1,1,'ArrayValued',true);


%yita
% Gt00t=int(tpy1*py1,y,-1,1);
% Gt11t=int(diff(tpy1,y)*diff(py1,y),y,-1,1);
% Gt02t=int(tpy1*diff(py1,y,2),y,-1,1);
% Gt20t=int(diff(tpy1,y,2)*py1,y,-1,1);
% Gt22t=int(diff(tpy1,y,2)*diff(py1,y,2),y,-1,1);

Gt00t=integral(matlabFunction(tpy1*py1),-1,1,'ArrayValued',true);
Gt11t=integral(matlabFunction(diff(tpy1,y)*diff(py1,y)),-1,1,'ArrayValued',true);
Gt02t=integral(matlabFunction(tpy1*diff(py1,y,2)),-1,1,'ArrayValued',true);
Gt20t=integral(matlabFunction(diff(tpy1,y,2)*py1),-1,1,'ArrayValued',true);
Gt22t=integral(matlabFunction(diff(tpy1,y,2)*diff(py1,y,2)),-1,1,'ArrayValued',true);


%梁
R00=tp00*p00;
R11=diff(tp00,z)*diff(p00,z);
%梁kesi
for i = 1:l
    if l==0
        break
    else
        F00{i}=eval(subs(R00,z,xl(1,i)));        %need to use numerical expression
        F11{i}=eval(subs(R11,z,xl(1,i)));
    end
end
%梁yita
for i=1:k;
    if k==0
        break
    else
        G00{i}=eval(subs(R00,z,yk(1,i)));
        G11{i}=eval(subs(R11,z,yk(1,i)));
    end
end

%板
%Ft22t*Gt00t

% Q1=[];
% for i=1:n;
%     Ai=[];
%     for j=1:n;
%         Bi=Ft22t(i,j)*ones(n,n);
%         A0=Bi.*Gt00t;
%         Ai=[Ai,A0];
%     end
%     Q1=[Q1;Ai];
% end

Q1=kron(Ft22t,Gt00t);       %this is simpler


%Ft00t*Gt22t

% Q2=[];
% for i=1:n;
%     Ai=[];
%     for j=1:n;
%         Bi=Ft00t(i,j)*ones(n,n);
%         A0=Bi.*Gt22t;
%         Ai=[Ai,A0];
%     end
%     Q2=[Q2;Ai];
% end

Q2=kron(Ft00t,Gt22t);


%Ft02t*Gt20t
% Q3=[];
% for i=1:n;
%     Ai=[];
%     for j=1:n;
%         Bi=Ft02t(i,j)*ones(n,n);
%         A0=Bi.*Gt20t;
%         Ai=[Ai,A0];
%     end
%     Q3=[Q3;Ai];
% end

Q3=kron(Ft02t,Gt20t);


%Ft20t*Gt02t
% Q4=[];
% for i=1:n;
%     Ai=[];
%     for j=1:n;
%         Bi=Ft20t(i,j)*ones(n,n);
%         A0=Bi.*Gt02t;
%         Ai=[Ai,A0];
%     end
%     Q4=[Q4;Ai];
% end
Q4=kron(Ft20t,Gt02t);

%Ft11t*Gt11t
% Q5=[];
% for i=1:n;
%     Ai=[];
%     for j=1:n;
%         Bi=Ft11t(i,j)*ones(n,n);
%         A0=Bi.*Gt11t;
%         Ai=[Ai,A0];
%     end
%     Q5=[Q5;Ai];
% end
Q5=kron(Ft11t,Gt11t);

%Ft00t*Gt00t
% Q6=[];
% for i=1:n;
%     Ai=[];
%     for j=1:n;
%         Bi=Ft00t(i,j)*ones(n,n);
%         A0=Bi.*Gt00t;
%         Ai=[Ai,A0];
%     end
%     Q6=[Q6;Ai];
% end
Q6=kron(Ft00t,Gt00t);


%梁
%Ft22t*G00
Q7=zeros(n^2,n^2);
% for i=1:k;
%     if k==0
%         break
%     else
%         A7=[];
%         for s=1:n;
%             Ai=[];
%             for j=1:n;
%                 Bi=Ft22t(s,j)*ones(n,n);
%                 A0=Bi.*G00{i};             
%                 Ai=[Ai,A0];
%             end
%             A7=[A7;Ai];
%         end
%         Q7=Q7+A7;
%     end
% end

for i=1:k
    if k==0
        break
    else
A7=kron(Ft22t,G00{i});    % faster 
        Q7=Q7+A7;
    end
end


%Ft11t*G11
% Q8=zeros(n^2,n^2);
% for i=1:k;
%     if k==0
%         break
%     else
%         A8=[];
%         for s=1:n;
%             Ai=[];
%             for j=1:n;
%                 Bi=Ft11t(s,j)*ones(n,n);
%                 A0=Bi.*G11{i};
%                 Ai=[Ai,A0];
%             end
%             A8=[A8;Ai];
%         end
%         Q8=Q8+A8;
%     end
% end
Q8=zeros(n^2,n^2);
for i=1:k
    if k==0
        break
    else
      A8=kron(Ft11t,G11{i});    % faster 
        Q8=Q8+A8;
    end
end

%F00*Gt22t
% Q9=zeros(n^2,n^2);
% for i=1:l;
%     if l==0
%         break
%     else
%         A9=[];
%         for s=1:n;
%             Ai=[];
%             for j=1:n;
%                 Bi=F00{i}(s,j)*ones(n,n);
%                 A0=Bi.*Gt22t;
%                 Ai=[Ai,A0];
%             end
%             A9=[A9;Ai];
%         end
%         Q9=Q9+A9;    
%     end
% end

Q9=zeros(n^2,n^2);
for i=1:l
    if l==0
        break
    else
     A9=kron(F00{i},Gt22t);    % faster 
        Q9=Q9+A9;    
    end
end

%F11*Gt11t
% Q10=zeros(n^2,n^2);
% for i=1:l;
%     if l==0
%         break
%     else
%         A10=[];
%         for s=1:n;
%             Ai=[];
%             for j=1:n;
%                 Bi=F11{i}(s,j)*ones(n,n);
%                 A0=Bi.*Gt11t;
%                 Ai=[Ai,A0];
%             end
%             A10=[A10;Ai];
%         end
%         Q10=Q10+A10;
%     end
% end

Q10=zeros(n^2,n^2);
for i=1:l
    if l==0
        break
    else
       A10=kron(F11{i},Gt11t);    % faster 
        Q10=Q10+A10;
    end
end

%Ft00t*G00
% Q11=zeros(n^2,n^2);
% for i=1:k;
%     if k==0
%         break
%     else
%         A11=[];
%         for s=1:n;
%             Ai=[];
%             for j=1:n;
%                 Bi=Ft00t(s,j)*ones(n,n);
%                 A0=Bi.*G00{i};
%                 Ai=[Ai,A0];
%             end
%             A11=[A11;Ai];
%         end
%         Q11=Q11+A11;
%     end
% end

Q11=zeros(n^2,n^2);
for i=1:k;
    if k==0
        break
    else
       A11=kron(Ft00t,G00{i});    % faster 
        Q11=Q11+A11;
    end
end

%F00*Gt00t
% Q12=zeros(n^2,n^2);
% for i=1:l;
%     if l==0
%         break
%     else
%         A12=[];
%         for s=1:n;
%             Ai=[];
%             for j=1:n;
%                 Bi=F00{i}(s,j)*ones(n,n);
%                 A0=Bi.*Gt00t;
%                 Ai=[Ai,A0];
%             end
%             A12=[A12;Ai];
%         end
%         Q12=Q12+A12;
%     end
% end
Q12=zeros(n^2,n^2);
for i=1:l;
    if l==0
        break
    else
       A12=kron(F00{i},Gt00t);    % faster 
        Q12=Q12+A12;
    end
end
%% stress matrix  sigma_xx
%stress matrix/vector of plate(middle point)---xx
x_loc=0.5;y_loc=0.5;
S_plt_mid_matrix1=-E*(-h/2)/(1-v^2)*((2/a)^2*eval(subs(diff(tp1,x,2),x,x_loc)*subs(diff(py1,y,0),y,y_loc))-...
    v*(2/b)^2*eval(subs(diff(tp1,x,0),x,x_loc)*subs(diff(py1,y,2),y,y_loc)));
S_plt_mid_vec1=reshape(S_plt_mid_matrix1',n*n,1);
% stress matrix/vector of beam(bottom of end point)--xx
x_loc=1;y_loc=max(yk);
S_beam_matrix1=-E*-(h/2+Hx+t2x)/(1-v^2)*((2/a)^2*eval(subs(diff(tp1,x,2),x,x_loc)*subs(diff(py1,y,0),y,y_loc)));
S_beam_vec1=reshape(S_beam_matrix1',n*n,1);

% %stress matrix/vector of plate---yy
% x_loc=0.5;y_loc=0.5;
% S_plt_mid_matrix2=-E*h/2/(1-v^2)*((2/b)^2*eval(subs(diff(tp1,x,0),x,x_loc)*subs(diff(py1,y,2),y,y_loc))+...
%     v*(2/a)^2*eval(subs(diff(tp1,x,2),x,x_loc)*subs(diff(py1,y,0),y,y_loc)));
% S_plt_mid_vec2=reshape(S_plt_mid_matrix2',n*n,1);
% % stress matrix/vector of beam--yy
% x_loc=xl;y_loc=1;
% S_beam_matrix2=-E*(h/2+Hy+t2y)/(1-v^2)*(2/b)^2*eval(subs(diff(tp1,x,0),x,x_loc)*subs(diff(py1,y,2),y,y_loc));
% S_beam_vec2=reshape(S_beam_matrix2',n*n,1);

%% displacement matrix middle of beam
x_loc=0.5;y_loc=max(yk);
D_beam_matrix1=eval(subs(diff(tp1,x,0),x,x_loc)*subs(diff(py1,y,0),y,y_loc));
D_beam_vec1=reshape(D_beam_matrix1',n*n,1);



%% stiffness and mass matrix
K=Q1+(r^4)*Q2+r^2*v*(Q3+Q4)+2*r^2*(1-v)*Q5+2/b/D*(E*Ix*Q7+r^2*G*Jx*Q8+r^3*E*Iy*Q9+r*G*Jy*Q10);
M=Q6+2/b/h*(Ax*Q11+1/r*Ay*Q12);

mod_K=4*b*D/(a^3);        %modification of K method 2
mod_M=a*b*P*h/4 ;         %modification of M method 2
A = K*mod_K;
B = M*mod_M;
% A = eval(K);
% B = eval(M);
[V,o]=eig(A,B);
% [V,o]=eig(A);
% w=(o*16*D/P/h/(a^4)).^(1/2)/2/pi;
w=sqrt(o)/(2*pi);           %frequency in Hz
[omega,w_order]=sort(diag(w));
mode=V(:,w_order);


figure
plot(1:9,omega(1:9),'o-');
xlabel('阶数');
ylabel('频率 (Hz)');
grid on;
omega(1:9)
t=toc


 %%  calculate mass and modification coefficient of K and M
%  mass=eval(sum(sum(M)));% total mass
 
 dt_default=0.1/max(omega);% default time step                                                                                                                                                                                     
 
 16*D/P/h/(a^4);     %modification of K method 1
 mod_K=4*b*D/(a^3);        %modification of K method 2
 mod_M=a*b*P*h/4;         %modification of M method 2
                                                   
% a*b*P*h/8*mass     % modified mass

%% plot mode shape
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

W_base=tp1*py1;
X_num=eval(subs(tp1,x,linspace(-a/2,a/2,31)));  % rows corresponds to i=1,2,...
Y_num=eval(subs(tpy1,y,linspace(-b/2,b/2,31))); % rows corresponds to j=1,2,...

W_base_num=cell(n,n);                       % basic function in numerical
for i=1:n
    for j=1:n
          W_base_num{i,j}=Y_num(j,:)'*X_num(i,:);       % can do this since X Y is independant
    end
end



%% plot vibration mode shape
    num=6;
    [X,Y] = meshgrid(linspace(-a/2,a/2,31),linspace(-b/2,b/2,31));
for k=1:num

Z_disp_vect=mode(:,k);      % displacement in 广义坐标
% Z_disp_vect=K_num\[1;zeros(224,1)]
    Z_disp{k}=zeros(31,31);
    for i=1:n
        for j=1:n
            Z_disp{k}=Z_disp{k}+W_base_num{i,j}*Z_disp_vect((i-1)*n+j,1);%displacement of whole plate 位移实际坐标
            
        end
    end

figure
    surf(X,Y,Z_disp{k});
    xlabel('X(m)');
    ylabel('Y(m)');
    zlabel('Z坐标(m)');
%     axis([-a/2,a/2,-b/2,b/2]);
%     zlim([0,1.5*Pm]);
    title(['Mode ',num2str(k)],['Frequency =',num2str(omega(k)),'Hz']);
end


%% save data
K_num=K;
M_num=M;
save (['data','.mat']);%保存所有数据
save (['K_M_matrix','.mat'],'K_num','M_num','n','mod_K','mod_M','W_base','W_base_num','S_plt_mid_vec1','S_beam_vec1','D_beam_vec1');%保存刚度矩阵、质量矩阵、阶数






