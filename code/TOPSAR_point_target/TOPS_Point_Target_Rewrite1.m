clear; clc;close all;
%%参数设置 
%%------------------------------------------------------------
%假设方位预处理使用的是补零升采样，方位后处理使用BAS方法
%%------------------------------------------------------------
c=physconst('LightSpeed');%光速
f0 = 9.65e9;              % 雷达工作频率GHz.
lambda=c/f0;              % 波长
Vr=7200;                  % 卫星速度
Fr = 20e6;                % 距离采样率
Fa = 3475;                % 见表4.1，PRF.

R_nc = 600e3;             %景中心斜距
Tb=0.48;                  %burst时间
omega_r=(3.415*pi)/180;   %波束转动速度，这里转成弧度
beta_bw =(0.33*pi)/180;     %方位向波束宽度，根据论文表4.1给出.
Lbeam = R_nc*tan(beta_bw);   %波束地面驻留长度，其中应该是用到了小角度近似
Vb=Vr+omega_r*R_nc;       %地面波束移动速度
Tsar = Lbeam/Vb;          %驻留时间
Td=Tsar;

Ka=-2*Vr^2/lambda/R_nc;  
Ba=abs(Ka*Tsar);            %点目标的多普勒带宽
Fa=3475;                    %方位采样率
PRF=Fa;
PRT=1/PRF;

Naz=ceil(Tb*Fa);            %方位向采样点数
t_slow=(-Naz/2:1:Naz/2-1)/Fa;   %慢时间

Tr = 10e-6;             
Kr = 20e12;
Br = Kr*Tr;         %距离向带宽200MHzs
Fr = 1.4*Br;        %距离向采样率
Y0=1000;            %当平台位于慢时间0时刻时，场景最短斜距为R_nc-Y0，最长斜距为R_nc+Y0
Rmin=R_nc-Y0;
Rmax=sqrt((R_nc+Y0)^2+(Tb*(Vb-Vr)/2)^2);
target_position=[0 R_nc+300;
                  0 R_nc;
                  2000 R_nc;
                  3000 R_nc;
                  1500 R_nc+250 ];%目标的位置
K=size(target_position,1);

Nrg=ceil(((Rmax-Rmin)/c*2+Tr)*Fr);
t_fast=(-Nrg/2:Nrg/2-1)/Fr+2*R_nc/c;%距离向快时间

echo_all=zeros(Naz,Nrg);
for k = 1:5
    A0=1;
    R=sqrt((t_slow*Vr-target_position(k,1)).^2+target_position(k,2)^2);%慢时间瞬时斜距，行向量
    tau=2*R/c;%时延
    wr=abs(ones(Naz,1)*t_fast-tau.'*ones(1,Nrg))<Tr/2;
    wa=((abs(t_slow*Vb-target_position(k,1))<Lbeam/2)'*ones(1,Nrg));
    echo=A0*exp(-j*4*pi/lambda*R.'*ones(1,Nrg)).*exp(1j*pi*Kr*(ones(Naz,1)*t_fast-tau.'*ones(1,Nrg)).^2).*wr.*wa;
    echo_all=echo_all+echo;%将k个目标加起来
end
figure;imagesc(real(echo_all));title("回波二维时域");
figure;imagesc(abs(ftx(echo_all)));title("回波RD域");
%% 方位变标预处理
% 解斜
t_slow=t_slow.';%转置
r_rot=-Vr/omega_r;      % 旋转中心到雷达平台飞行航迹的最短斜距,按照论文中的r坐标轴的取值来定义
k_rot=-2*Vr^2/(lambda*r_rot);%是正值
echo_all=echo_all.*exp(-1j*pi*k_rot*t_slow.^2*ones(1,Nrg));
figure;imagesc(real(echo_all));title("方位解斜后二维时域");
% 频域补零
ka=-2*Vr^2/(lambda*R_nc);%方位向调频率，论文公式4.16，r取景中心斜距R_nc
Bd=abs(Ka*Td);%点目标有效多普勒带宽
Bf=2*Vr/lambda*tan(beta_bw);%方位波束带宽
Bs=Tb*abs(k_rot)+Bf;

PRF_NEW = ceil(Bs);
N_az=ceil(PRF_NEW*(PRF/k_rot/2))*2;%新的方位向采样点数，PRF/Krot是方位PRF时间，ceil内/2再在外面*2是为了保证N_az是偶数
%?这里用的是PRF时间，但是，实际的Burst时间Tb要比PRF时间大很多，所以应该不对
N_az=ceil(Tb*PRF_NEW/2)*2;
sig=zeros(N_az,Nrg);
echo_all=ftx(echo_all);%方位向变到多普勒域

sig(N_az/2-Naz/2:N_az/2+Naz/2-1,:)=echo_all;%补零
t_slow=(-N_az/2:N_az/2-1)/PRF_NEW;%新的慢时间
t_slow=t_slow.';
f_dc=0;
fa=(-N_az/2:N_az/2-1)/N_az*PRF_NEW+f_dc;
fa=fa.';

echo_all=sig;
figure;imagesc(real(echo_all));title("频域补零后的RD数据");
figure;imagesc(real(iftx(echo_all)));title("频域补零后的时域数据");
% 恢复多普勒历程
echo_all=iftx(echo_all);
echo_all=echo_all.*exp(1j*pi*k_rot*t_slow.^2*ones(1,Nrg));
figure;imagesc(real(echo_all));title("方位恢复多普勒历程后二维时域");
%% CS设置
r_ref=R_nc;%参考距离
R=t_fast/2*c;
beta_fa=sqrt(1-(lambda*fa/(2*Vr)).^2);%式子4.23
a_fa=1./beta_fa-1;%CS因子，见式子4.24
k_fa_r=1./(1/(-Kr)-2*lambda*((beta_fa.^2-1)./(c^2*beta_fa.^3))*R);%4.22，RD域的调频率,与fa,r有关
k_fa_ref=1./(1/(-Kr)-(2*lambda*R_nc*(beta_fa.^2-1))./(c^2*beta_fa.^3));%式子里的r都换成景中心斜距R_nc
R_fa_r=(1+a_fa)*R;%随着fa变化的斜距，式子4.21，r是斜距，那么理解成不同距离门的回波？
R_fa_ref=R_nc*(1+a_fa);%随着fa变化的斜距，式子4.21
%% 距离向chirp scaling
echo_all=ftx(echo_all);
echo_all=echo_all.*exp(-1j*pi*k_fa_ref.*a_fa*ones(1,Nrg).*(ones(N_az,1)*t_fast-2*R_fa_ref/c*ones(1,Nrg)).^2);
figure;imagesc(real(echo_all));title("距离向CS后RD域");
%% 距离脉压和一致RCMC
echo_all=fty(echo_all);
fr = ( -Nrg/2 : Nrg/2-1 )*( Fr/Nrg);
echo_all=echo_all.*exp(-1j*pi*1./(k_fa_ref.*(1+a_fa))*fr.^2).*exp(1j*4*pi*r_ref/c*a_fa*fr);
echo_all=ifty(echo_all);
figure;imagesc(abs(echo_all));title("距离向脉压和一致RCMC后RD域");
%% H3 相位补偿
ascl_fa=(1+a_fa).*sqrt(1-(lambda*f_dc/(2*Vr)).^2);
echo_all=echo_all.*exp(1j*4*pi*k_fa_ref.*ascl_fa.*(1+a_fa).^2/c^2./(1+ascl_fa)*(R-r_ref).^2);
figure;imagesc(abs(echo_all));title("H3后RD域");
%%
r0=R_nc;%参考距离
r_scl_r=R_nc/r_rot*(r_rot-R)/(1-(R_nc/r_rot));%见Prats的论文
echo_all=echo_all.*exp(1j*4*pi/lambda*(beta_fa-1)*R).*exp(1j*pi*lambda/(2*Vr^2)*(fa).^2*r_scl_r); 
k_rot1=-2*Vr^2./(lambda*(r_rot-R)/(1-R_nc/r_rot));
echo_all=iftx(echo_all);
echo_all=echo_all.*exp(-1j*pi*t_slow.^2*k_rot1);
echo_all=ftx(echo_all);
ka2=-2*Vr^2./(lambda*r_scl_r)-k_rot1;
echo_all=echo_all.*exp(1j*pi*(fa).^2*(1./ka2));
echo_all=iftx(echo_all);
figure;imagesc(abs(echo_all));
target_1 = target_analysis( echo_all,Fr,Fa,Vb);%注意时间的选取
%%---------------------------------------------------------
%以下内容是按照陈保坤论文和Prats论文中的公式写的，其中，尝试使用了解斜之后，方位采样之前降采样处理，但是还是结果不对
% %% H4
% r_scl0=r_rot*(1-PRF/PRF_NEW);
% r_rot_r=(r_rot-R)./(1-r_scl0/r_rot);
% r_scl=r_scl0/r_rot*r_rot_r;
% k_scl=-2*Vr^2/lambda./r_scl;
% echo_all=echo_all.*exp(1j*4*pi/lambda*(beta_fa-1)*R).*exp(1j*pi*fa.^2*(1./k_scl));
% figure;imagesc(abs(echo_all));title("H4后RD域");
% %% H5
% echo_all=iftx(echo_all);
% k_rot=-2*Vr^2/lambda./r_rot_r;
% echo_all=echo_all.*exp(-1j*pi*t_slow.^2*k_rot);
% figure;imagesc(abs(echo_all));title("H5后时域");
% %% 降采样
% echo_all=ftx(echo_all);
% Ns=ceil(Bf/PRF_NEW*N_az/2)*2;
% sig=echo_all(N_az/2-Ns/2+1:N_az/2+Ns/2,:);
% echo_all=sig;
% t_slow=(-Ns/2:Ns/2-1)/Bf;%新的慢时间
% t_slow=t_slow.';
% fa=(-Ns/2:Ns/2-1)/Ns*Bf+f_dc;
% fa=fa.';
% 
% %% H6
% Keff=k_scl-k_rot;
% echo_all=ftx(echo_all);
% echo_all=echo_all.*exp(1j*pi*fa.^2*(1./Keff)).*(rectpuls(fa,1900)*ones(1,Nrg));
% echo_all=iftx(echo_all);
% figure;imagesc(abs(echo_all));title("H6后时域");
%%---------------------------------------------------------------------

%%---------------------------------------------------------------
%以下公式是按照TOPS-Mode Raw Data Processing Using Chirp Scaling
%Algorithm的内容写的，但是不知道为嘛结果不对，存疑
% %% H3,4
% alpha=0.9;
% alpha_a=alpha*R_nc./R;%行向量
% ascl_fa=(1+a_fa).*sqrt(1-(lambda*f_dc/(2*Vr)).^2);
% echo_all=echo_all.*exp(1j*4*pi/lambda*(beta_fa-1)*R).*exp(1j*pi*lambda*r_rot/2/Vr^2*fa.^2*alpha_a)...
%     .*exp(1j*4*pi*k_fa_ref.*ascl_fa.*(1+a_fa).^2/c^2./(1+ascl_fa)*(R-r_ref).^2);
% figure;imagesc(abs(echo_all));title("H3,4后RD域");
% %% H4,4 de-rotation
% echo_all=iftx(echo_all);
% echo_all=echo_all.*exp(-1j*pi*k_rot*t_slow.^2*R_nc*(1./R));
% figure;imagesc(abs(echo_all));title("H4,4后二维时域");
% %% H5,4
% echo_all=ftx(echo_all);
% ka_new=2*Vr^2/lambda*1./(R+r_rot-alpha_a*r_rot);
% echo_all=echo_all.*(exp(1j*pi*fa.^2*1./(ka_new-alpha_a*k_rot)).*(rectpuls(fa,12000)*ones(1,Nrg)));%方位压缩
% figure;imagesc(abs(echo_all));title("H5,4后RD域");
% %% H6,4
% echo_all=iftx(echo_all);
% echo_all=echo_all.*exp(1j*pi*t_slow.^2*1./(ka_new-alpha_a*k_rot)*alpha^2);
% figure;imagesc(abs(echo_all));title("H6,4后二维时域");
%%---------------------------------------------------------------------