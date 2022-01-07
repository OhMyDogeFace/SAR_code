%% 尝试子孔径分块处理方法
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
Naz = 1657;                 % 距离线数（即数据矩阵，行数）——这里修改为1024。不要设置得太大，否则方位向会出现问题
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
%% 设定仿真点的位置
%先用一个burst进行处理，我们仍然使用斜距-方位距离坐标系来定义点目标的位置
%定义五个距离向位置相同，方位向位置不同的点目标先，以景中心的目标两边对称
target_position=[R_nc 0;
                R_nc 60;
                R_nc-50 0;
                R_nc 80;
                R_nc 80;
                            ];%距离向12m可以分辨出来,方位向要到32m可以分辨出来
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
for k =1:2
    wa=(abs(t_slow) <= Tb/2);% 先对burt时间内的数据，加窗
    wa1=(abs(t_slow - nc_target(k)) <= Td);% 先不加sinc窗，仅仅加矩形窗;
    theta=atan(Vr*(t_slow-nc_target(k))/target_position(k,1))+omega_r*t_slow;% 式4.31，用于求回波包络
    % 注意这里多了一个omega_r*t_slow,主要是考虑到方位向的波束转动，从而导致不同的方位向天线图加权
    % 通过这里我们可以明白，theta代表的是理想天线方向图中，波束中心与视线的夹角
    wa1=(sinc(0.886*theta/beta_bw)).^2.*(abs(t_slow - nc_target(k)) <= Td);
    for i =1:Naz
        A0=1;
        R_n=sqrt(target_position(k,1)^2+(target_position(k,2)-t_slow(i)*Vr)^2);%计算瞬时的斜距
        wr=rectpuls(t_fast-2*R_n/c,Tr);%矩形窗包络,经验证，合理正确
        echo(i,:)=A0*wa(i)*wa1(i)*wr.*exp(-1j*4*pi*f0*R_n/c).*exp(-1j*pi*Kr*(t_fast-2*R_n/c).^2);%一个回波信号
    end
    echo_all=echo_all+echo;%将k个目标加起来     
end
figure(1);
imagesc(1:Nrg,1:Naz,real(echo_all));title('原始回波二维时域实部');
%% 分块处理
% 每个数据块的长度要大于波束在每个点的驻留时间，且彼此之间重叠5%~10%，时域上划分成数据块
T_sub=0.3*Td;%子数据块的时间?子数据块的时间该如何取,感觉要保证每个数据块的多普勒频率不大于PRF=3475Hz
%为了方便仿真，分成三个子孔径，这是预先设置好的
N_sub=floor(T_sub*Fa);%子数据块的数据点数
N_overlap=7;%
N=ceil(Naz/N_sub);
%见论文Spotlight SAR Data Processing Using the Frequency Scaling
%Algorithm的式子(39)和(40)关于子孔径频率的计算，TOPSAR中要将atan部分换成omega_r*ta
sub_echo=zeros(N_sub,Nrg,N);%三维矩阵
sub_t_slow=zeros(N,N_sub);
sub_f_nc=zeros(N,3);%N×3的矩阵，第一列为左多普勒频率，第二列为右多普勒频率，第三列为多普勒中心频率
for i =1 : N
    
sub_echo(:,:,i)=echo_all((1+(i-1)*(N_sub-N_overlap)):(N_sub+(i-1)*(N_sub-N_overlap)),:);%将原始仿真数据划分为N个部分，其中有重叠N_overlap个数据点
sub_t_slow(i,:)=t_slow((1+(i-1)*(N_sub-N_overlap)):(N_sub+(i-1)*(N_sub-N_overlap)));%同样将慢时间轴也这样划分成N个子孔径
sub_f_nc(i,1)=2*Vr*sin(omega_r*sub_t_slow(i,1)-beta_bw/2)/lambda;%计算每个子孔径的左右多普勒频率
sub_f_nc(i,2)=2*Vr*sin(omega_r*sub_t_slow(i,end)+beta_bw/2)/lambda;
sub_f_nc(i,3)=(sub_f_nc(i,1)+sub_f_nc(i,2))/2;%每个子孔径的多普勒频率
end
%后续操作和原先一样
%% 结合方位变标的CS算法
r_ref=R_nc;%参考距离，先写成景中心斜距，目前假设的斜视角为0
R=t_fast/2*c;%距离向时间轴上面，不同的距离门内，不同的斜距
sub_echo_RD1=zeros(size(sub_echo));
sub_echo_RD2=zeros(size(sub_echo));
sub_echo_RD3=zeros(size(sub_echo));
sub_echo_RD4=zeros(size(sub_echo));
for m =1:N
    fa=fftshift((-N_sub/2:N_sub/2-1)/N_sub*(sub_f_nc(m,2)-sub_f_nc(m,1))+sub_f_nc(m,3));  %方位向频率轴的选取，多普勒中心频率为各自的
    beta_fa=sqrt(1-(lambda*fa/(2*Vr)).^2);%式子4.23
    a_fa=1./beta_fa-1;%CS因子，见式子4.24
    k_fa_r=1./(1/Kr-2*lambda*((beta_fa.^2-1)./(c^2*beta_fa.^3)).'*R);%4.22，RD域的调频率,与fa,r有关
    k_fa_ref=1./(1/Kr-(2*lambda*R_nc*(beta_fa.^2-1))./(c^2*beta_fa.^3));%式子里的r都换成景中心斜距R_nc
    R_fa_r=(1+a_fa).'*R;%随着fa变化的斜距，式子4.21，r是斜距，那么理解成不同距离门的回波？
    R_fa_ref=R_nc*(1+a_fa);%随着fa变化的斜距，式子4.21
    %信号变到RD域
    sub_echo_RD=fft(sub_echo(:,:,m),N_sub,1);%方位向进行fft,零频在两端，此时频率轴与fa是对应起来的
    sub_echo_RD1(:,:,m)=zeros(size(sub_echo_RD));
    % 距离向chirp scaling操作
    for n = 1:N_sub
        H1=exp(-1j*pi*k_fa_ref(n)*a_fa(n)*(t_fast-2*R_fa_ref(n)/c).^2);%式子4.25,实际上公式中的r是参考距离r_ref，但是这里简化都写成R_nc，即景中心斜距
                                                               %因为先写的斜视角为0,k_fa_应该也是k_fa_ref
        sub_echo_RD1(n,:,m)=sub_echo_RD(n,:).*H1;
    end
%     figure;
%     imagesc(1:Nrg,1:N_sub,real(fftshift(sub_echo_RD,1)));title('距离向chir scaling操作前RD域实部');
%     figure;
%     imagesc(1:Nrg,1:N_sub,real(fftshift(sub_echo_RD1(:,:,m),1)));title('距离向chir scaling操作后RD域实部');
    % 距离向脉压和一致RCMC
    sub_echo_2f=fft(sub_echo_RD1(:,:,m),Nrg,2);%距离向fft，换到二维频域，此时距离向方位向频率轴的零频都在两端
    sub_echo_2f1=zeros(size(sub_echo_2f));
    fr= (-Nrg/2 : Nrg/2-1 )*( Fr/Nrg );%距离向频率轴,为了将距离频率轴对应起来，
    t_ref=( -Nrg/2 : (Nrg/2-1) )/Fr; 
    for n = 1:N_sub
        H2=exp(-1j*pi*fr.^2/(k_fa_ref(n)*(1+a_fa(n)))).*exp(j*4*pi*r_ref/c*a_fa(n)*fr);
        H2=fftshift(H2);
        sub_echo_2f1(i,:)=sub_echo_2f(i,:).*H2;
    end
    sub_echo_RD2(:,:,m)=ifft(sub_echo_2f1,Nrg,2);%距离向ifft,变到RD域
     figure;
     imagesc(1:Nrg,1:N_sub,abs(fftshift(sub_echo_RD2(:,:,m),1)));title('距离向脉压和一致RCMC后RD域实部');
    % 相位校正
    % 相位校正函数，是随着fa变化的，同时需要引入不同距离们（即不同斜距代表的r）与r_ref之间的差delta_phi_fa=-4*pi/(c)^2*(k_fa_ref.').*(a_fa.').*(1+a_fa).'*(R-r_ref*ones(size(t_fast))).^2;%写成矩阵的形式
    for n = 1:Nrg
        delta_phi_fa=4*pi/(c)^2*k_fa_ref.*a_fa.*(1+a_fa)*(R(i)-r_ref).^2;%写成矩阵的形式
        H3=exp(1j*delta_phi_fa);%相位校正函数式子4.30
        sub_echo_RD3(:,n,m)=sub_echo_RD2(:,n,m).*H3.';
    end
    % 方位变标,文章中说H5的引进是为了移除双曲线型相位并引入一个二次线性调频相位,就是方位变标
    f_sdc=0;%场景多普勒中心，假设为0
    r0=R_nc;%参考距离
    r_scl_r=(2*R+omega_r*R.^2/Vr-r0)./(1+omega_r*R./Vr);%式子4.38
    r_scl_r=R_nc/r_rot*(r_rot-R)/(1-(R_nc/r_rot));%见Prats的论文
    for n = 1:Nrg
        H4=exp(1j*4*pi/lambda*R(n)*(beta_fa-1)).*exp(1j*pi*lambda*r_scl_r(n)./(2*Vr^2).*(fa).^2); 
        sub_echo_RD4(:,n,m)=sub_echo_RD3(:,n,m).*H4.';
    end
end
%% 子孔径拼接
% %接下来的子孔径拼接，对于重叠部分的处理，我觉得应该要去除，因为后面还有解斜操作，如果不去除就直接拼接，方位向时间会有冗余部分
 echo_all2=zeros(size(echo_all));
 sub_echo1=zeros(size(sub_echo));
for m = 1:N-1
    sub_echo1(:,:,m)=ifft(sub_echo_RD4(:,:,m),N_sub,1);%每个子孔径数据方位ifft
%      figure;
%      imagesc(abs(sub_echo1(:,:,m)));
    echo_all2(1+(m-1)*(N_sub-N_overlap):m*(N_sub-N_overlap),:)=sub_echo1(1:N_sub-N_overlap,:,m);
end
sub_echo1(:,:,N)=ifft(sub_echo_RD4(:,:,N),N_sub,1);
echo_all2(1+(N-1)*(N_sub-N_overlap):end,:)=sub_echo1(1:N_sub,:,m);
%存疑，感觉方位拼接得有点问题，
%% 方位时域去斜
k_rot1=2*Vr^2./(lambda*(Vr/omega_r+R-r_scl_r));%新的多普勒中心变化率，见式子4.40
k_rot1=-2*Vr^2./(lambda*(r_rot-R)/(1-R_nc/r_rot));
echo_all3=zeros(size(echo_all2));
for i =1:Nrg
    hd1=exp(-1j*pi*k_rot1(i)*t_slow.^2);
    echo_all3(:,i)=echo_all2(:,i).*hd1.';
end
echo_all_RD1=fft(echo_all3,Naz,1);%方位向fft，变到RD域
%% 方位压缩
ka2=-2*Vr^2./(lambda*r_scl_r)-k_rot1;
echo_all_RD2=zeros(size(echo_all_RD1));
fa=fftshift((-Naz/2:Naz/2-1)/Naz*Fa);  %方位向频率轴的选取，多普勒中心频率为0
for i =1:Nrg
    H6=exp(1j*pi*(fa-f_sdc).^2/ka2(i)).*fftshift(rectpuls(fftshift(fa),400));%方位匹配滤波函数,注意这里仍然要取符号，不是正号
    echo_all_RD2(:,i)=echo_all_RD1(:,i).*H6.';
end
echo_image=ifft(echo_all_RD2,Naz,1);
%%
r_rot_r=(r_rot-R)/(1-R_nc/r_rot);
K_t_r=-2*Vr^2./(lambda*(r_rot_r-r_scl_r));
echo_image1=zeros(size(echo_image));
for i = 1:Nrg
    H7=exp(1j*pi*K_t_r(i)*(1-R_nc/r_rot)^2*t_slow.^2);
   echo_image1(:,i)=echo_image(:,i).*H7.';
end
figure;
imagesc(1:Nrg,1:Naz,real(fftshift(echo_all_RD2,1)));title('方位向压缩后RD域实部');
figure;
imagesc(1:Nrg,1:Naz,abs(echo_image1));title('点目标');
figure;
mesh(1:Nrg,1:Naz,abs(echo_image1));title('三维图');

NN=48;
target_1 = target_analysis( echo_image1,Fr,Fa,Vr);
%% 先子孔径成像？
% k_rot1=2*Vr^2./(lambda*(Vr/omega_r+R-r_scl_r));%新的多普勒中心变化率，见式子4.40
% k_rot1=-2*Vr^2./(lambda*(r_rot-R)/(1-R_nc/r_rot));
% sub_echo2=zeros(size(sub_echo));
% sub_echoRD5=zeros(size(sub_echo));
% sub_echoRD6=zeros(size(sub_echo));
% sub_echoimage=zeros(size(sub_echo));
% sub_echoimage1=zeros(size(sub_echo));
% for m = 1:N
%     
%     for n = 1:Nrg
%         hd1=exp(-1j*pi*k_rot1(n)*sub_t_slow(m,:).^2);
%         sub_echo2(:,n,m)=sub_echo1(:,n,m).*hd1.';
%     end
%     ka2=-2*Vr^2./(lambda*r_scl_r)-k_rot1;
%     fa=fftshift((-N_sub/2:N_sub/2-1)/N_sub*(sub_f_nc(m,2)-sub_f_nc(m,1))+sub_f_nc(m,3));  %方位向频率轴的选取，多普勒中心频率为各自的
%     sub_echoRD5(:,:,m)=fft(sub_echo2(:,:,m),N_sub,1);
%     for n =1:Nrg
%         H6=exp(1j*pi*(fa-sub_f_nc(m,3)).^2/ka2(n));%方位匹配滤波函数,注意这里仍然要取符号，不是正号
%         sub_echoRD6(:,n,m)=sub_echoRD5(:,n,m).*H6.';
%     end
%     sub_echoimage(:,:,m)=ifft(sub_echoRD6(:,:,m),N_sub,1);
%     r_rot_r=(r_rot-R)/(1-R_nc/r_rot);
%     K_t_r=-2*Vr^2./(lambda*(r_rot_r-r_scl_r));
%     for n = 1:Nrg
%         H7=exp(1j*pi*K_t_r(n)*(1-R_nc/r_rot)^2*sub_t_slow(m,:).^2);
%         sub_echoimage1(:,n,m)=sub_echoimage(:,n,m).*H7.';
%     end
%     figure;
%     mesh(abs(sub_echoimage1(:,:,m)));
% end
