clear; clc;close all;
%%参数设置 
%%------------------------------------------------------------------
%方位预处理：两步式成像
%方位后处理：重调频，方位匹配滤波
%--------------------------------------------------------------------
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
target_position=[0 R_nc;
                  3000 R_nc;
                  -3000 R_nc;
                  -300 R_nc;
                  -6000 R_nc+500 ];%目标的位置
K=size(target_position,1);

Nrg=ceil(((Rmax-Rmin)/c*2+Tr)*Fr);
t_fast=(-Nrg/2:Nrg/2-1)/Fr+2*R_nc/c;%距离向快时间

echo_all=zeros(Naz,Nrg);
for k = 3
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
%% TOPS模式因子的设置
Rr=c*t_fast/2;%随着快时间变化的的不同距离门,行向量
A=1+omega_r*Rr/Vr;%滑动因子，见星载TOPSAR模式研究论文中的4.1
Rrot=Rr/(1-A);%旋转中心到航迹的距离.随着距离门变化而变化?可直接取论文中的,是带入化简的结果
Rrot=-Vr/omega_r;
Krot = -2*Vr^2/lambda/Rrot; %%多普勒中心频率变化率

Bf=abs(2*Vr/lambda*tan(beta_bw));  %方位波束带宽
Bb=Tb*Krot+Bf;              %Burst时间总带宽
Mpre = Bb/PRF;              %频域扩展比例

t_slow=t_slow.';
deramp = exp(1j*pi*(-Krot)*t_slow.^2)*ones(1,Nrg);
echo_all=echo_all.*deramp;
figure;imagesc(real(echo_all));title("方位去斜后的二维时域");
figure;imagesc(real(ftx(echo_all)));title("方位去斜后的RD域");

PRF_NEW = Mpre*PRF;     %新的PRF
N_az=ceil(PRF_NEW*(PRF/Krot/2))*2;%新的方位向采样点数，PRF/Krot是方位PRF时间，ceil内/2再在外面*2是为了保证N_az是偶数
sig = zeros(N_az,Nrg);
% echo_all=ftx(echo_all);%方位向fft!!!!!!!!!!!!!!!!!!!!

sig(N_az/2-Naz/2:N_az/2+Naz/2-1,:)=echo_all;%补零
t_slow=(-N_az/2:N_az/2-1)/PRF_NEW;%新的慢时间
t_slow=t_slow.';
f_dc=0;
fa=(-N_az/2:N_az/2-1)/N_az*PRF_NEW+f_dc;
fa=fa.';

echo_all=sig;
figure;imagesc(real(echo_all));title("时域补零后的RD数据");
figure;imagesc(real(ftx(echo_all)));title("时域补零后的RD域域数据");

%% 方位向两步处理
echo_all=ftx(echo_all);
echo_all=echo_all.*(exp(1j*pi*(-Krot)*t_slow.^2)*ones(1,Nrg));
figure,imagesc(real(echo_all));title('相位补偿');

echo_all=ftx(echo_all);
echo_all=echo_all.*(exp(1i*pi*fa.^2/(-Krot))*ones(1,Nrg));%恢复多普勒历程
figure,imagesc(real(echo_all));title('恢复多普勒历程');
%% CS设置
R_ref=R_nc;%参考目标的斜距
fn_ref=f_dc;
D_fn_Vr = sqrt(1-lambda^2.*(fa).^2./(4*Vr^2));% 大斜视角下的徙动因子，式7.17,列向量
D_fn_ref_Vr = sqrt(1-lambda^2*fn_ref^2/(4*Vr^2));% 参考频率fn_ref处的徙动因子，是常数。
K_src = 2*Vr^2*f0^3*D_fn_Vr.^3./(c*R_ref*(fa).^2);% 式6.22,二次距离压缩滤波器的调频率，列向量
                                                  % 此处R取R_ref，见7.33下方的描述
Km = Kr./(1-Kr./K_src);% 这是变换到距离多普勒域的距离调频率，列向量
                       % 使用 R_ref 处的值
%% 距离向chirpscaling                       
echo_all=echo_all.*(exp(1j*pi*Km.*(D_fn_ref_Vr./D_fn_Vr-1)*(t_fast-2*R_ref./(c.*D_fn_ref_Vr)).^2));
figure,imagesc(real(echo_all));title('CS相位因子相乘');
%% 一致距离徙动校正和距离压缩SRC
fr = ( -Nrg/2 : Nrg/2-1 )*( Fr/Nrg);
echo_all=fty(echo_all);
echo_all=echo_all.*exp(1j*pi*D_fn_Vr./(Km*D_fn_ref_Vr)*fr.^2).*exp(1j*4*pi/c*(1./D_fn_Vr-1/D_fn_ref_Vr)*R_ref*fr);

echo_all=ifty(echo_all);
figure;imagesc(abs(echo_all));title("距离校正后RD域");

%% 方位压缩和相位校正
phase_azimuth=exp(1j*4*pi*R_nc*f0*D_fn_Vr/c).*ones(1,Nrg)...%方位向匹配滤波器
    .*(exp(-1j*4*pi.*Km/(c^2).*(1-D_fn_Vr/D_fn_ref_Vr).*(R_nc./D_fn_Vr-R_ref./D_fn_Vr).^2)*ones(1,Nrg));%附加相位校正
echo_all=echo_all.*phase_azimuth;
figure;imagesc(abs(echo_all));title("方位压缩和相位校正RD域");

%%
K_scl=Krot*Ka./(Ka-Krot);%见《星载TOPSAR模式研究》式4.29，根据几何关系可以得出
echo_all=echo_all.*exp(1j*pi*fa.^2*ones(1,Nrg)./K_scl);%这里的fa也是混叠的,这应该是频域去斜函数,变成直的了
echo_all=iftx(echo_all);%变到时域
figure;imagesc(abs(echo_all));%如果没有论文中的G2补偿相位项，不知为什么会出现方位多个目标，也许是解斜函数频率混叠？
echo_all=echo_all.*exp(1j*pi*t_slow.^2*ones(1,Nrg).*K_scl);%利用SPECAN原理，单频信号FFT后会变成sinc信号
echo_all=ftx(echo_all);%
figure;imagesc(abs(echo_all));
 
target_1 = target_analysis( echo_all,Fr,PRF,Vb);

