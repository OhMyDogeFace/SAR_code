%% 尝试TOPSAR点目标仿真，参照RDA的代码
%方位预处理：解斜频域补零来解决方位混叠
%方位后处理：方位变标CS
clear; close all;
%% 参数设置
%根据自己的设置来进行变动
R_nc = 600e3;               % 景中心斜距,目标被波束中心照射时的雷达据目标的斜距,看作是场景中心到平台的距离.
Vr = 7200;                  % 速度m/s.
Tr = 10e-6;                % 发射脉冲时宽us
Kr = 1.5e12;                 % 距离调频率Hz/s
f0 = 9.65e9;                % 雷达工作频率GHz.
% BW_dop = 80;                % 多普勒带宽Hz，也可以根据方位向调频率乘以合成孔径时间算出来=79.8508
Fr = 20e6;                  % 距离采样率.
Fa = 3475;                  % 见表4.1，PRF.
Naz = 1024;                 % 距离线数（即数据矩阵，行数）——这里修改为1024。不要设置得太大，否则方位向会出现问题
Nrg = 512;                  % 距离线采样点数（即数据矩阵，列数）这里设置的是采样点数。
% sita_r_c = (21.9*pi)/180;	% 波束斜视角，3.5 度，这里转换为弧度.
% 假设是正侧视
c = physconst('LightSpeed');% 光速.
Tb = 0.48;                  % burst长度0.48s.
Td = 0.1131;                % 驻留时间.可以理解为一个合成孔径时间？
BW_range = Kr*Tr;           % 距离向带宽
lambda = c/f0;               % 波长.
% fnc = 2*Vr*sin(sita_r_c)/lamda;     % 多普勒中心频率，根据公式（4.33）计算。

beta_bw =(0.33*pi)/180;     %方位向波束宽度，根据论文表4.1给出.
La_real=0.886*lambda/beta_bw;
La = 0.886*R_nc*lambda/La_real;% 合成孔径长度Ls(4.49)Vs/Vg约等于1.

% Mamb = round(fnc/Fa);       % 多普勒模糊
omega_r=(3.415*pi)/180;     %波束转动速度，这里转成弧度
r_rot=-Vr/omega_r;      % 旋转中心到雷达平台飞行航迹的最短斜距,按照论文中的r坐标轴的取值来定义
Vf=omega_r*R_nc+Vr;     % 波束在地面移动的速度
Lbeam = R_nc*tan(beta_bw);   %波束地面驻留长度，其中应该是用到了小角度近似
%% 设定仿真点的位置
%先用一个burst进行处理，我们仍然使用斜距-方位距离坐标系来定义点目标的位置
%定义五个距离向位置相同，方位向位置不同的点目标先，以景中心的目标两边对称
target_position=[R_nc 0;
                R_nc 33;
                R_nc+600 -4000;
                R_nc-400 5000;
                R_nc+400 -6000;
                R_nc+300 -2400;
                R_nc+250 2400;
                R_nc -8000];%距离向12m可以分辨出来,方位向要到32m可以分辨出来
nc_target=zeros(size(target_position));
for i = 1:size(nc_target,1)
    nc_target(i)=atan(target_position(i,2)/(target_position(i,1)-r_rot))/omega_r;%我是根据目标的旋转角大小值来计算波束中心穿越时刻的
end
%% 定义时间轴
t_fast=(-Nrg/2:Nrg/2-1)/Fr+2*target_position(1,1)/c;%距离时间轴
t_slow=(-Naz/2:Naz/2-1)/Fa;%慢时间轴，飞机飞行的时间轴
sita_c=omega_r*t_slow;      %波束中心与零多普勒线之间的夹角，波束中心斜视角，从负到正,实际上应该根据burst持续时间Tb定义，但是为了方便，后面再重新加窗
%% 一个子带，一个burst的回波生成
echo = zeros(Naz,Nrg);
echo_all = zeros(Naz,Nrg);
for k = 1:1
    A0=1;
    R_n=sqrt(target_position(k,1)^2+(target_position(k,2)*ones(1,Naz)-t_slow*Vr).^2);%每个方位向上的时间点到目标的距离,行向量
    tau_mtx=2*(R_n/c).'*ones(1,Nrg);%每个慢时间点的时延,矩阵形式
    wr=abs(ones(Naz,1)*t_fast-tau_mtx)<Tr/2;%距离窗
    theta=atan(Vr*(t_slow-nc_target(k))/target_position(k,1))*pi/180+omega_r*t_slow;% 式4.31，用于求回波包络,弧度制
    wa1=((sinc(0.886*theta/beta_bw)).^2.*(abs(t_slow - nc_target(k)) <= Td)).'*ones(1,Nrg);%方位窗
    wa=((abs(t_slow*Vf-target_position(k,2))<Lbeam/2)'*ones(1,Nrg));
    echo=A0*(exp(-1j*4*pi*f0*R_n/c)).'*ones(1,Nrg).*exp(-1j*pi*Kr*(ones(Naz,1)*t_fast-tau_mtx).^2).*wa.*wr;
    echo_all=echo_all+echo;%将k个目标加起来     
end
figure(1);
imagesc(1:Nrg,1:Naz,abs(echo_all));title('原始回波二维时域实部');
%% 方位变标预处理
k_rot=-2*Vr^2/(lambda*r_rot);%是正值
hd=exp(-1j*pi*k_rot*t_slow.^2).*rectpuls(t_slow,Tb);%去斜函数，这里是混叠的，见卡明书中的图9.3
echo_all_1=zeros(size(echo_all));
for i = 1:Nrg
   echo_all_1(:,i)=echo_all(:,i).*hd.'; 
end
figure(2);
imagesc(1:Nrg,1:Naz,real(echo_all_1));title('回波方位向乘以去斜函数后二维时域实部');
%再求Bs—方位向Burst信号总带宽 Bs=Tb*k_rot+Bf;式子4.13
%Bf-方位向瞬时带宽（方位波束带宽）Bf=A_r*Bd; %式子4.12
%Bd-点目标有效多普勒带宽Bd=ka*Td;%根据图4.3TOPSAR回波方位时频，或者就是驻留时间乘方位调频率
sita_0=0;%Burst中心时刻方位向波束指向
ka=-2*Vr^2*(cos(sita_0))^3/(lambda*R_nc);%方位向调频率，论文公式4.16，r取景中心斜距R_nc
A_r=(-r_rot+R_nc)/-r_rot;%公式4.1，方位分辨率改变因子
Bd=abs(ka*Td);
Bf=A_r*Bd;
Bs=Tb*k_rot+Bf;
echo_all_1_F=fft(echo_all_1,Naz,1);%方位向变换到频域，为补零做准备，此时方位向零频在两端
echo_all_1_F=fftshift(echo_all_1_F,1);%fftshift，此时零频在中间
N_az=ceil(Bs/Fa*Naz);%补零后，新的方位向数据点数
echo_all_F_zeros_padding=[zeros((N_az-Naz)/2,Nrg);echo_all_1_F;zeros((N_az-Naz)/2,Nrg)];%方位向补零，零频在中间
echo_all_F_zeros_padding=fftshift(echo_all_F_zeros_padding,1);%零频在两端
echo_all_2=ifft(echo_all_F_zeros_padding,N_az,1);%逆fft，此时方位向数据点数为N_az
figure(3);
imagesc(1:Nrg,1:N_az,real(echo_all_2));title('回波方位向乘以去斜函数、频域补零后二维时域实部');
echo_all_3=zeros(size(echo_all_2));
t_slow_2=(-N_az/2:N_az/2-1)/Bs;%此时采样率是Bs
%用去斜函数的共轭还原多普勒历程
hd1=exp(1j*pi*k_rot*t_slow_2.^2).*rectpuls(t_slow_2,Tb);
for i =1:Nrg
    echo_all_3(:,i)=echo_all_2(:,i).*hd1.';%.*(abs(t_slow_2 - nc_target(1)) <= Td).'; 升采样之后的信号还需要限定驻留时间Td吗？
                                           %应该不用，因为方位向上错落着不同方位距离的目标回波，不需要也不能够进行限定驻留时间Td
end
figure(4);
imagesc(1:Nrg,1:N_az,real(echo_all_3));title('回波方位向乘以去斜函数、频域补零、方位向乘去斜共轭后二维时域实部');
%% 结合方位变标的CS算法
r_ref=R_nc;%参考距离，先写成景中心斜距，目前假设的斜视角为0
R=t_fast/2*c;%距离向时间轴上面，不同的距离门内，不同的斜距
fa=fftshift((-N_az/2:N_az/2-1)/N_az*Bs);  %方位向频率轴的选取，多普勒中心频率为0
beta_fa=sqrt(1-(lambda*fa/(2*Vr)).^2);%式子4.23
a_fa=1./beta_fa-1;%CS因子，见式子4.24
k_fa_r=1./(1/Kr-2*lambda*((beta_fa.^2-1)./(c^2*beta_fa.^3)).'*R);%4.22，RD域的调频率,与fa,r有关
k_fa_ref=1./(1/Kr-(2*lambda*R_nc*(beta_fa.^2-1))./(c^2*beta_fa.^3));%式子里的r都换成景中心斜距R_nc
R_fa_r=(1+a_fa).'*R;%随着fa变化的斜距，式子4.21，r是斜距，那么理解成不同距离门的回波？
R_fa_ref=R_nc*(1+a_fa);%随着fa变化的斜距，式子4.21
%先将信号变到RD域
echo_all_RD=fft(echo_all_3,N_az,1);%方位向进行fft,零频在两端，此时频率轴与fa是对应起来的
echo_all_RD1=zeros(size(echo_all_RD));
%% 距离向chirp scaling操作
for i = 1:N_az
    H1=exp(-1j*pi*k_fa_ref(i)*a_fa(i)*(t_fast-2*R_fa_ref(i)/c).^2);%式子4.25,实际上公式中的r是参考距离r_ref，但是这里简化都写成R_nc，即景中心斜距
                                                               %因为先写的斜视角为0,k_fa_应该也是k_fa_ref
    echo_all_RD1(i,:)=echo_all_RD(i,:).*H1;
end
figure(5);
imagesc(1:Nrg,1:N_az,real(fftshift(echo_all_RD1,1)));title('距离向chir scaling操作后RD域实部');
%% 距离向脉压和一致RCMC
echo_all_2f=fft(echo_all_RD1,Nrg,2);%距离向fft，换到二维频域，此时距离向方位向频率轴的零频都在两端
echo_all_2f1=zeros(size(echo_all_2f));
fr= (-Nrg/2 : Nrg/2-1 )*( Fr/Nrg );%距离向频率轴,为了将距离频率轴对应起来，
t_ref=( -Nrg/2 : (Nrg/2-1) )/Fr; 
for i = 1:N_az
    H2=exp(-1j*pi*fr.^2/(k_fa_ref(i)*(1+a_fa(i)))).*exp(j*4*pi*r_ref/c*a_fa(i)*fr);
    %注意这里论文的一个错误，流程图里面没有标注共轭
    H2=fftshift(H2);
    echo_all_2f1(i,:)=echo_all_2f(i,:).*H2;
end
echo_all_RD2=ifft(echo_all_2f1,Nrg,2);%距离向ifft,变到RD域
figure(6);
imagesc(1:Nrg,1:N_az,abs(fftshift(echo_all_RD2,1)));title('距离向脉压和一致RCMC后RD域实部');
%% 相位校正
% 相位校正函数，是随着fa变化的，同时需要引入不同距离们（即不同斜距代表的r）与r_ref之间的差delta_phi_fa=-4*pi/(c)^2*(k_fa_ref.').*(a_fa.').*(1+a_fa).'*(R-r_ref*ones(size(t_fast))).^2;%写成矩阵的形式
echo_all_RD3=zeros(size(echo_all_RD2));
for i = 1:Nrg
    delta_phi_fa=4*pi/(c)^2*k_fa_ref.*a_fa.*(1+a_fa)*(R(i)-r_ref).^2;%写成矩阵的形式
    H3=exp(1j*delta_phi_fa);%相位校正函数式子4.30
    echo_all_RD3(:,i)=echo_all_RD2(:,i).*H3.';
end
figure(7);
imagesc(1:Nrg,1:N_az,real(fftshift(echo_all_RD3,1)));title('相位校正后RD域实部');
%% 天线幅度校正
 f_sdc=0;%场景多普勒中心，假设为0
%    H4=(sinc(1*(fa-f_sdc).*beta_fa/(2*Vr))).^(-2);%天线幅度校正函数
%    echo_all_RD4=zeros(size(echo_all_RD3));
%     for i =1:Nrg
%         echo_all_RD4(:,i)=echo_all_RD3(:,i).*H4.';
%     end
%     figure(8);
%    imagesc(1:Nrg,1:N_az,real(fftshift(echo_all_RD4,1)));title('天线幅度校正后RD域实部');
%% 方位变标,文章中说H5的引进是为了移除双曲线型相位并引入一个二次线性调频相位,就是方位变标
%见A SAR processing algorithm for TOPS imaging mode一文中关于azimuth scaling的描述
r0=R_nc;%参考距离
r_scl_r=(2*R+omega_r*R.^2/Vr-r0)./(1+omega_r*R./Vr);%式子4.38
r_scl_r=R_nc/r_rot*(r_rot-R)/(1-(R_nc/r_rot));%见Prats的论文
echo_all_RD5=zeros(size(echo_all_RD3));
for i =1:Nrg
   H5=exp(1j*4*pi/lambda*R(i)*(beta_fa-1)).*exp(1j*pi*lambda*r_scl_r(i)./(2*Vr^2).*(fa).^2); 
   echo_all_RD5(:,i)=echo_all_RD3(:,i).*H5.';
end
figure(9);
imagesc(1:Nrg,1:N_az,real(fftshift(echo_all_RD5,1)));title('方位变标后RD域实部');
%% 方位时域去斜
k_rot1=2*Vr^2./(lambda*(Vr/omega_r+R-r_scl_r));%新的多普勒中心变化率，见式子4.40
k_rot1=-2*Vr^2./(lambda*(r_rot-R)/(1-R_nc/r_rot));
echo_all_4=ifft(echo_all_RD5,N_az,1);%方位向ifft，变到二维时域
echo_all_5=zeros(size(echo_all_4));
for i =1:Nrg
    hd1=exp(-1j*pi*k_rot1(i)*t_slow_2.^2);
    echo_all_5(:,i)=echo_all_4(:,i).*hd1.';
end
echo_all_RD6=fft(echo_all_5,N_az,1);%方位向fft，变到RD域
%% 进行降采样Bw>=Bd1=Ba=Bf
echo_all_RD6=fftshift(echo_all_RD6,1);
echo_all_RD7=echo_all_RD6((N_az/2-round(Bf/Bs*N_az)/2+1):(N_az/2+round(Bf/Bs*N_az/2)),:);
echo_all_RD7=fftshift(echo_all_RD7,1);
figure(10);
imagesc(1:Nrg,1:N_az,real(fftshift(echo_all_RD7,1)));title('方位时域去斜后、降采样RD域实部');
N_az_new=round(Bf/Bs*N_az);%此时采样率是Bf
fa1=fftshift(-N_az_new/2:N_az_new/2-1)/N_az_new*Bf;  %方位向频率轴的选取，多普勒中心频率为0
%% 方位压缩
ka2=-2*Vr^2./(lambda*r_scl_r)-k_rot1;
echo_all_RD8=zeros(size(echo_all_RD7));
for i =1:Nrg
    H6=exp(1j*pi*(fa1-f_sdc).^2/ka2(i)).*fftshift(rectpuls(fftshift(fa1),1900));%方位匹配滤波函数,注意这里仍然要取符号，不是正号
    echo_all_RD8(:,i)=echo_all_RD7(:,i).*H6.';
end
echo_image=ifft(echo_all_RD8,N_az_new,1);
% figure(11);
% imagesc(1:Nrg,1:N_az_new,real(fftshift(echo_all_RD8,1)));title('方位向压缩后RD域实部');
% figure(12);
% imagesc(1:Nrg,1:N_az_new,abs(echo_image));title('点目标');
% figure(13);
% mesh(1:Nrg,1:N_az_new,abs(echo_image));title('三维图');
%%
t_slow_3=(-N_az_new/2:N_az_new/2-1)/Bf;
r_rot_r=(r_rot-R)/(1-R_nc/r_rot);
K_t_r=-2*Vr^2./(lambda*(r_rot_r-r_scl_r));
echo_image1=zeros(size(echo_image));
for i = 1:Nrg
    H7=exp(1j*pi*K_t_r(i)*(1-R_nc/r_rot)^2*t_slow_3.^2);
   echo_image1(:,i)=echo_image(:,i).*H7.';
end

figure(11);
imagesc(1:Nrg,1:N_az_new,real(fftshift(echo_all_RD8,1)));title('方位向压缩后RD域实部');
figure(12);
imagesc(1:Nrg,1:N_az_new,abs(echo_image1));title('点目标');
figure(13);
mesh(1:Nrg,1:N_az_new,abs(echo_image1));title('三维图');

%% 点目标升采样分析

NN=48;
target_1 = target_analysis( echo_image1,Fr,Fa,Vf);