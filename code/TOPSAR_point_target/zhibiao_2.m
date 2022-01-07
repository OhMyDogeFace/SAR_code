function [PSLR,ISLR,IRW] = zhibiao_2(x,s_number,T)
%
% ��Ժ��� zhibiao(x,s_number,T) ���иĽ����Ľ���������ҪΪ��
% �ڼ���3dB������Ӧ������ʱ������������ �ֱ��ʣ�IRW��
% ���� zhibiao() �в��õ����ٽ�ȡ���İ취���ⲻ׼ȷ��
% �������õķ����ǽ��� 3dB �����������������Բ�ֵ�����õ���׼ȷ��3dB������Ӧ�����ꡣ
%
% ����������ź�x���������� s_number��T���źŵ�ʱ�򳤶ȡ�
% �ú���������� x �ķ�ֵ�԰��(PSLR)�������԰��(ISLR),����ֱ��ʣ�IRW����
soo = x;                % �ź�x
N_buling = s_number;    % ��������

soo_abs = abs(soo);     % soo��ģ
[C,I] = max(soo_abs);   % ����� soo��ģ �е����ֵ C��λ�� I��
y = soo_abs.^2;         % �����ƽ���� y = soo^2��

x1 = 0;
while (soo_abs(I-x1-1)-soo_abs(I-x1))<0 %�����ߵ�ֵ���ұߵ�С
    M1 = x1;%������ƫ�Ƶĵ�Ԫ������M1�У�ֱ����������һ�԰�ʱ����ʱ������ߵ���͵���I-x1-1
    x1 = x1+1;
end
x2 = 0;
while (soo_abs(I+x2+1)-soo_abs(I+x2))<0%
    M2 = x2;%ͬ��I+M2+1���ұ���͵�
    x2 = x2+1;
end

P1 = I-1-M1;            % ������԰�ֽ�㣬��ߵ������� P1��
P2 = I+1+M2;            % ������԰�ֽ�㣬�ұߵ������ P2��

[D_left,Q_left] = max(soo_abs(1,1:P1));     % ����԰ֵ꣬Ϊ D_left��λ��Ϊ Q_left������ߵ���һ������
[D_right,Q_right] = max(soo_abs(1,P2:end)); % ����԰ֵ꣬Ϊ D_right��λ��Ϊ Q_right�����ұߵ���һ������
D = max(D_left,D_right);    % �Ƚ���ߺ��ұ������е����ֵ���õ������԰��������԰ֵ꣬Ϊ D��

PSLR = 20*log10(D/C);                       % ��ֵ�԰��
ISLR = 10*log10((sum(y(1,1:P1))+sum(y(1,P2:end)))/sum(y(1,P1:P2)));% �����԰�ȣ�10lg((P_total-P_main)/P_main))

%%%%%%%%%%%%%%%%%%%%%%%  �������� IRW  %%%%%%%%%%%%%%%%%%%%%%%%%

M = ( 10^(-3/20) )*C;       % 3dB �����ĺ���ȡֵ,M�ķ�ֵ��C��0.7079�������ʾ��Ƕ���֮1��
% ������Ϊ�����ҳ���ú���ֵ��ӽ���ֵ�Ĵ�С�����ꡣ
[value_left,index_left]=min(abs(soo_abs(1:I)-M));%�ҳ��ӷ�ֵ�����M��ӽ��ĵ�����ꣻ
[value_right,index_right]=min(abs(soo_abs(I+1:end)-M));%�ҳ��ӷ�ֵ�����M��ӽ��ĵ�����ꣻ
width=index_right-index_left+I;
c = physconst("LightSpeed");% ���� c=3e8 m/s��
IRW = T/N_buling*width*c/2;% ע����SAR�У��ֱ����� C*T/2,����T�������ȡ�
% IRW_realΪͼ��ֱ��ʣ�ԭ����width�ĵ�λ�ǲ��������һ��ͼ��ֱ��ʵ�λȡ m��Ҫת����
% ע�⵽���������õ��� N_buling����ΪƵ������ЧΪ����������������ߣ���������Ӧ�ø���Ϊ N_buling��
