function [PSLR,ISLR,IRW] = zhibiao_2(x,s_number,T)
%
% 针对函数 zhibiao(x,s_number,T) 进行改进，改进的内容主要为：
% 在计算3dB点所对应的坐标时——用来计算 分辨率（IRW）
% 函数 zhibiao() 中采用的是临近取整的办法，这不准确。
% 下面利用的方法是将离 3dB 最近的两个点进行线性插值，来得到更准确的3dB点所对应的坐标。
%
% 输入变量：信号x，采样点数 s_number，T是信号的时域长度。
% 该函数用来求解 x 的峰值旁瓣比(PSLR)，积分旁瓣比(ISLR),距离分辨率（IRW）。
soo = x;                % 信号x
N_buling = s_number;    % 采样点数

soo_abs = abs(soo);     % soo的模
[C,I] = max(soo_abs);   % 求输出 soo的模 中的最大值 C，位置 I；
y = soo_abs.^2;         % 输出的平方， y = soo^2。

x1 = 0;
while (soo_abs(I-x1-1)-soo_abs(I-x1))<0 %如果左边的值比右边的小
    M1 = x1;%将向左偏移的单元数记在M1中，直到滑动到第一旁瓣时，此时主瓣左边的最低点是I-x1-1
    x1 = x1+1;
end
x2 = 0;
while (soo_abs(I+x2+1)-soo_abs(I+x2))<0%
    M2 = x2;%同理，I+M2+1是右边最低点
    x2 = x2+1;
end

P1 = I-1-M1;            % 主瓣和旁瓣分界点，左边的坐标是 P1。
P2 = I+1+M2;            % 主瓣和旁瓣分界点，右边的左边是 P2。

[D_left,Q_left] = max(soo_abs(1,1:P1));     % 最大旁瓣，值为 D_left，位置为 Q_left。（左边的那一个）。
[D_right,Q_right] = max(soo_abs(1,P2:end)); % 最大旁瓣，值为 D_right，位置为 Q_right。（右边的那一个）。
D = max(D_left,D_right);    % 比较左边和右边两者中的最大值，得到两侧旁瓣中最大的旁瓣，值为 D。

PSLR = 20*log10(D/C);                       % 峰值旁瓣比
ISLR = 10*log10((sum(y(1,1:P1))+sum(y(1,P2:end)))/sum(y(1,P1:P2)));% 积分旁瓣比，10lg((P_total-P_main)/P_main))

%%%%%%%%%%%%%%%%%%%%%%%  以下是求 IRW  %%%%%%%%%%%%%%%%%%%%%%%%%

M = ( 10^(-3/20) )*C;       % 3dB 带宽处的函数取值,M的幅值是C的0.7079倍，功率就是二分之1。
% 下面是为了求找出与该函数值最接近的值的大小和坐标。
[value_left,index_left]=min(abs(soo_abs(1:I)-M));%找出从峰值左边与M最接近的点和坐标；
[value_right,index_right]=min(abs(soo_abs(I+1:end)-M));%找出从峰值左边与M最接近的点和坐标；
width=index_right-index_left+I;
c = physconst("LightSpeed");% 光速 c=3e8 m/s。
IRW = T/N_buling*width*c/2;% 注意在SAR中，分辨率是 C*T/2,其中T是脉冲宽度。
% IRW_real为图像分辨率，原来的width的单位是采样间隔。一般图像分辨率单位取 m，要转换。
% 注意到采样点数用的是 N_buling，因为频域补零后等效为升采样，采样率提高，采样点数应该跟正为 N_buling。
