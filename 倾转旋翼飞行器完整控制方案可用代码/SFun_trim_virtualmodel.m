%���ܣ�����S_Functionģ�飬�Զ���������ת�����������ѧϵͳ��ģ��
function [sys,x0,str,ts] = SFun_trim_virtualmodel(t,x,u,flag)
xinit = [0;0;0;0;0];
switch flag,
case 0,
   [sys,x0,str,ts] = mdlInitializeSizes(xinit);                            % ��ʼ������
case 1,
    sys = mdlDerivatives(t,x,u);                                 % ����״̬��������ģ��
case 3,                             
    sys = mdlOutputs(t,x);                                                 % ���ģ��
case { 2, 4, 9 }                       
    sys = [];                                                              % δ�����־
otherwise
    error(['Unhandled flag = ',num2str(flag)]);                            % �������
end

 
function [sys,x0,str,ts] = mdlInitializeSizes(xinit)
sizes = simsizes;           % ϵͳĬ�����ã�����δ��ʼ������
sizes.NumContStates = 5;    % ����״̬��������
sizes.NumDiscStates = 0;    % ��ɢ״̬����0
sizes.NumOutputs = 5;       % ����������������ϵͳ�����״̬��
sizes.NumInputs = 2;        % �������������������������������ţ�
sizes.DirFeedthrough = 0;   % ֱͨ���� 0 ��ʾϵͳ�����ɱ�������ڲ����������
sizes.NumSampleTimes = 1;   % �������ڸ���
sys = simsizes(sizes);      % ��size��ʼ������򴫸�sys
str = [];                   % �����ַ�������
ts = [0 0];                 % -1�̳������źŲ�������
x0 = xinit;                 % ״̬�����ĳ�ֵ


function sys = mdlDerivatives(t,x,u)
[Fxt,Fyt,M_z] = Force_Moment_Compute(t,x,u);                               %������������

global tao_online; tao = tao_online;
global t_online; t = t_online;

t_wing = 0; t_hover = 0;
global t_tran;                                                             %����ʱ��
global taov0;global v0_down;                                               %�ն�����     
v0_up = 23.5; v_tran = v0_up - v0_down;                                    %78������ 
%S3 Vx����
if t == 0 || t >= t_tran || tao > taov0-1  %״̬�����ĵ����������յ㲻�ɵ�����Ҫ���⴦��(ˮƽ�����ȼ����˶�)
    accelerate_vx = -1.4289; 
else
    accelerate_vx = -v_tran / t_tran * sin(pi * t / t_tran) / (sqrt(1-(2 * tao / taov0 - 1)^2));
end

p1 = 7.7191e-08; p2 = 6.8331e-06; p3 = -0.00018488; p4 = 0.0012743; p5 = -0.0034819; p6 = 0.02417;
accelerate_theta = 6 * p1 * t^5 + 5 * p2 * t^4 + 4 * p3 * t^3 + 3 * p4 * t^2 + 2 * p5 * t + p6; %�ڶ��׶���ƽ��ָ����theta������Ҫ��
%accelerate_theta = 0;                                                      %��һ�׶���ƽ����theta������Ҫ��
accelerate_vy = 0;
accelerate_height = 0;
accelerate_q = 0;

dx1 = (Fxt*cos(x(3)) - Fyt*sin(x(3)))/2 - accelerate_vx;                   %״̬���̣�����ϵͳ״̬΢��          
dx2 = (Fxt*sin(x(3)) + Fyt*cos(x(3)))/2 - accelerate_vy;
dx3 = x(4) - accelerate_theta;
dx4 = M_z/0.0514 - accelerate_q;
dx5 = x(2) - accelerate_height;
sys = [dx1;dx2;dx3;dx4;dx5];


function sys = mdlOutputs(t,x)
sys = x;


function [Fxt,Fyt,M_z] = Force_Moment_Compute(t,x,u)
global tao_online;                                                         %�̶��� ȫ��������ת�Ǹ���
global tao_trimed;
%tao = tao_trimed;%0;%tao_online;
tao = tao_online;
tau = tao * pi / 180;                                                      %����ת��                                               
cy_alfa = 0.0247; cy_wz = 2.91;   cy_deltaz = 0.00737;                     %���������������
mz_alfa = -2.04;  mz_wz = -5.865; mz_deltaz = -0.0184;
range = 1;%0.6;
cy_alfa = range*cy_alfa; cy_wz = range*cy_wz; cy_deltaz = range*cy_deltaz; %�����㶯
mz_alfa = range*mz_alfa; mz_wz = range*mz_wz; mz_deltaz = range*mz_deltaz;

velocity = sqrt((x(1))^2 + (x(2))^2); theta = x(3); alfa = theta - atan(x(2)/x(1)); %������ٺ�ӭ��
x(1) = velocity; x(2) = alfa;

C_y = cy_alfa*x(2)*180/pi + 0.1503 + cy_wz*0.21*x(4)/x(1) + cy_deltaz*u(1);%��������������
C_y1 = cy_alfa*x(2)*180/pi + 0.1503;
C_x = 0.3328*(C_y1)^2 - 0.0635*C_y1 + 0.0233;          
Q = (0.5*1.225*x(1)^2)*0.233*C_x;                                          
Y = (0.5*1.225*x(1)^2)*0.233*C_y;  

T = 14.75*u(2) - 0.819;                                                    %������������
Fxtj = -Q*cos(x(2)) + Y*sin(x(2)) - 2*9.8*sin(x(3));
Fxtx = cos(tau)*T*2;
Fytj = Q*sin(x(2)) + Y*cos(x(2)) - 2*9.8*cos(x(3));
Fytx = sin(tau)*T*2;

v_hl = -136.95*u(2)^4 + 385.2*u(2)^3 - 404.46*u(2)^2 + 194.49*u(2) + 39.741;%դ������������
alfa_s1 = atan(x(1)*sin(x(2) + tau)/(v_hl + x(1)*cos(x(2) + tau)));
alfa_s = alfa_s1 * 180 / pi;
v_s = sqrt(v_hl^2 + x(1)^2 + 2*v_hl*x(1)*cos(x(2) + tau));
C_ys = (0.0131*alfa_s + 0.0519)/20;
C_xs = 0.0294*C_ys + 0.0079/20;
Q_s = 0.5*1.225*0.233*0.5*v_s^2*C_xs;                                      
Y_s = 0.5*1.225*0.233*0.5*v_s^2*C_ys;
Fxts = (-cos(tau - alfa_s1)*Q_s + sin(alfa_s1 - tau)*Y_s)*2;
Fyts = (-sin(tau - alfa_s1)*Q_s + cos(tau - alfa_s1)*Y_s)*2;

Fxt = Fxtj + Fxts + Fxtx;                                                  %�ػ�������ϵ�������
Fyt = Fytj + Fyts + Fytx;

dalfa = x(4) - (Fxt*sin(x(2)) + Fyt*cos(x(2)))/(2*x(1));                   %ӭ��alfa�ĵ���
mz1 = -0.1624*C_y1 + 0.0395 + mz_wz*(0.21*x(4)/x(1)) + mz_alfa*(0.21*dalfa/x(1)) + mz_deltaz*u(1);

M_zj = (0.5*1.225*x(1)^2)*0.233*0.21*mz1;                                  
M_zx = -T*cos(tau)*0.055*2;
M_zs = -2*0.055*(-cos(tau - alfa_s1)*Q_s + sin(alfa_s1 - tau)*Y_s) + 2*0.0055*(-sin(tau - alfa_s1)*Q_s + cos(tau - alfa_s1)*Y_s) + (0.5*1.225*v_s^2)*0.233*0.21*(-0.032*C_ys + 0.013/20);
M_z = (M_zj + M_zx + M_zs);                                                %�ػ�������ϵ���������