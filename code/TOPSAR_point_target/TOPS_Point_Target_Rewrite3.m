clear; clc;close all;
%%参数设置 
%%------------------------------------------------------------
%假设方位预处理使用的是mosaic，方位后处理使用基于SPECAN的ECS方法
%但是不知为何，方位向距离变大时，聚焦效果不好
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
                  -2000 R_nc;
                  -6000 R_nc+500 ];%目标的位置
K=size(target_position,1);

Nrg=ceil(((Rmax-Rmin)/c*2+Tr)*Fr);
t_fast=(-Nrg/2:Nrg/2-1)/Fr+2*R_nc/c;%距离向快时间

echo_all=zeros(Naz,Nrg);
for k = 1:4
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
%% 方位mosaic处理
Bd=abs(Ka*Td);%点目标有效多普勒带宽
Bf=2*Vr/lambda*tan(beta_bw);%方位波束带宽
r_rot=-Vr/omega_r;      % 旋转中心到雷达平台飞行航迹的最短斜距,按照论文中的r坐标轴的取值来定义
k_rot=-2*Vr^2/(lambda*r_rot);%是正值
Bs=Tb*abs(k_rot)+Bf;

Mf=ceil(Bs/PRF);%需要扩展的份数
t_slow=t_slow.';%转置

echo_all=ftx(echo_all);
N_az_new=Mf*Naz;
sig=zeros(N_az_new,Nrg);
for i = 1:Mf
   sig(1+(i-1)*Naz:Naz+(i-1)*Naz,:)=echo_all; 
end
echo_all=sig;%频域上拓展副本
%% 频域去斜
fa=(-N_az_new/2:N_az_new/2-1)/N_az_new*Mf*PRF;%新的方位频域范围
fa=fa.';
echo_all=echo_all.*exp(1j*pi*fa.^2/k_rot);
%% 方位时域截取
echo_all=iftx(echo_all);
Rr=c*t_fast/2;%随着快时间变化的的不同距离门,行向量
A=1+omega_r*R_nc/Vr;%滑动因子，见星载TOPSAR模式研究论文中的4.1
Tw=1.35*Tsar*A/(A-1);%这里时间的截取，对分辨率和旁瓣模糊有很大影响
Ns=ceil(Tw/Tb*N_az_new/2)*2;
sig=echo_all(N_az_new/2-Ns/2+1:N_az_new/2+Ns/2,:);
echo_all=sig;%方位时域的截取 此时方位时间范围为[-Tw/2,Tw/2],方位采样点数为Ns,频域范围为[-Mf*PRF/2,Mf*PRF/2]
%方位时间范围变小，频域范围变大
figure;imagesc(abs(echo_all));title("mosaic预处理后的回波二维时域图")
%% CS设置
t_slow=(-Ns/2:Ns/2-1)/(Mf*PRF);%新的慢时间
t_slow=t_slow.';
f_dc=0;
fa=(-Ns/2:Ns/2-1)/Ns*Mf*PRF+f_dc;
fa=fa.';

%  Nb=ceil(Ns/2/1.1)*2;
%  sig=zeros(Ns,Nrg);
% sig(Nb/2-Ns/2+1:Nb/2+Ns/2,:)=echo_all;
% sig(Ns/2-Nb/2+1:Ns/2+Ns/2,:)=echo_all(Ns/2-Nb/2+1:Ns/2+Ns/2,:);
% echo_all=sig;
% t_slow=(-Nb/2:Nb/2-1)/(Mf*PRF);
% t_slow=t_slow.';
% fa=(-Nb/2:Nb/2-1)/Nb*Mf*PRF+f_dc;
% fa=fa.';

r_ref=R_nc;%参考距离
R=t_fast/2*c;
beta_fa=sqrt(1-(lambda*fa/(2*Vr)).^2);%式子4.23
a_fa=1./beta_fa-1;%CS因子，见式子4.24
k_fa_r=1./(1/(-Kr)-2*lambda*((beta_fa.^2-1)./(c^2*beta_fa.^3))*R);%4.22，RD域的调频率,与fa,r有关
k_fa_ref=1./(1/(-Kr)-(2*lambda*R_nc*(beta_fa.^2-1))./(c^2*beta_fa.^3));%式子里的r都换成景中心斜距R_nc
R_fa_r=(1+a_fa)*R;%随着fa变化的斜距，式子4.21，r是斜距，那么理解成不同距离门的回波？
R_fa_ref=R_nc*(1+a_fa);%随着fa变化的斜距，式子4.21
%% 距离向chirp scaling
echo_all=ftx(echo_all);sig=ones(Ns,Nrg);
%echo_all(1:Ns/2-1002,:)=zeros(Ns/2-1002,Nrg);echo_all(Ns/2+1002+1:Ns,:)=zeros(Ns/2-1002,Nrg);%加不加影响不大好像
echo_all=echo_all.*exp(-1j*pi*k_fa_ref.*a_fa*ones(1,Nrg).*(ones(Ns,1)*t_fast-2*R_fa_ref/c*ones(1,Nrg)).^2);
figure;imagesc(real(echo_all));title("距离向CS后RD域");
%% 距离脉压和一致RCMC
echo_all=fty(echo_all);
fr = ( -Nrg/2 : Nrg/2-1 )*( Fr/Nrg);
echo_all=echo_all.*exp(-1j*pi*1./(k_fa_ref.*(1+a_fa))*fr.^2).*exp(1j*4*pi*r_ref/c*a_fa*fr);
echo_all=ifty(echo_all);
figure;imagesc(abs(echo_all));title("距离向脉压和一致RCMC后RD域");
%% 基于SPECAN的ECS算法
delta_x=1/(Mf*PRF)*Vb;%方位向采样间隔。等于方位向采样时间*速度，我觉得应该用波束移动速度
r_scl=2*Vr^2*delta_x*Ns/lambda/PRF/Vb;%变标距离
r_scl=(2*R_nc+omega_r*Rr.^2/Vr-R_nc)./(1+omega_r*Rr/Vr);%行向量
ascl_fa=(1+a_fa).*sqrt(1-(lambda*f_dc/(2*Vr)).^2);
H32=exp(1j*4*pi/lambda*(beta_fa-1)*Rr).*exp(1j*pi*lambda/2/Vr^2*fa.^2*(r_rot+r_scl))...
    .*exp(1j*4*pi*k_fa_ref.*ascl_fa.*(1+a_fa).^2/c^2./(1+ascl_fa)*(Rr-r_ref).^2);
echo_all=echo_all.*H32;
figure;imagesc(abs(echo_all));
%% H4,2
k_scl=-2*Vr^2/lambda./r_scl;%行向量
echo_all=iftx(echo_all);
echo_all=echo_all.*exp(-1j*pi*t_slow.^2*k_scl);
figure;imagesc(abs(echo_all));
%% H5,2
echo_all=ftx(echo_all);
echo_all=echo_all.*exp(1j*pi*lambda/2/Vr^2*fa.^2*(Rr-r_rot));

figure;imagesc(abs(echo_all));
target_1 = target_analysis( echo_all,Fr,Fa,Vb);%注意时间的选取