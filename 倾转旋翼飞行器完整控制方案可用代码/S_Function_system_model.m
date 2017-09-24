%���ܣ�����S_Functionģ�飬�Զ���������ת�����������ѧϵͳ��ģ��
%{
��������S_Function_system_model
����˵����t ��ǰʱ��
         x ��ǰϵͳ��״̬�������� x1�ٶ�  x2ӭ��  x3������  x4����������  x5�߶�
         u ��ǰϵͳ�Ŀ��������
         flag ģ�ͽӿڲ�����ָʾ�ò�ִ�еĹ���
         xinit ��ʼֵ����ΪS_F�ĸ��Ӳ���
         tao_trim ��ƽ�����ת��
���ý����sys ��ӳϵͳģ����Ϣ������������ֵȡ����flag
         x0  ϵͳ��ʼ״̬
         str State ordering strings which is generally specified as [].
         ts  ��ӳ����ʱ����������ɢ
example:[23.5 0.0398 0.0398 0 20]' or [3 0.0265 0.0265 0 20]'
 %}
function [sys,x0,str,ts] = S_Function_system_model(t,x,u,flag,xinit)
   switch flag,
   case 0,
       [sys,x0,str,ts] = mdlInitializeSizes(xinit);                        % ��ʼ������
   case 1,
        sys = mdlDerivatives(t,x,u);                                       % ����״̬��������ģ��
   case 3,                             
        sys = mdlOutputs(t,x);                                             % ���ģ��
   case { 2, 4, 9 }                       
        sys = [];                                                          % δ�����־
   otherwise
        error(['Unhandled flag = ',num2str(flag)]);                        % �������
   end
  
   
%{
��������mdlInitializeSizes
���ܣ�����S_F��ʼ��������ϵͳ������Ϣ
����˵����xinit ϵͳ��ʼ״̬
example:
 %}
function [sys,x0,str,ts] = mdlInitializeSizes(xinit)
sizes = simsizes;           % ϵͳĬ�����ã�����δ��ʼ������
sizes.NumContStates = 5;    % ����״̬��������
sizes.NumDiscStates = 0;    % ��ɢ״̬����0
sizes.NumOutputs = 16;      % ����������������ϵͳ�����״̬��
sizes.NumInputs = 2;        % �������������������������������ţ�
sizes.DirFeedthrough = 0;   % ֱͨ���� 0 ��ʾϵͳ�����ɱ�������ڲ����������
sizes.NumSampleTimes = 1;   % �������ڸ���
sys = simsizes(sizes);      % ��size��ʼ������򴫸�sys
str = [];                   % �����ַ�������
ts = [0 0];                 % -1�̳������źŲ�������
x0 = xinit;                 % ״̬�����ĳ�ֵ


%{
��������mdlDerivatives 
���ܣ�����ϵͳ����״̬�����ĵ��� 
����˵����t ��ǰʱ�� x ��ǰ״̬ u ��ǰ�������� ���ý����sys
example:
%}
function sys = mdlDerivatives(t,x,u)
if t == 0                                                                  %ϵͳ��ʼ��
    %System_Global_init(); 
end

[tao,v_tao] = systemcorridor(t);                                           %�ƶ���������                                                         
global tao_online; tao_online = tao;                                                          
global v_tao_online; v_tao_online = v_tao;
global t_online; t_online = t;
 
model = 0;                                                                 %���м�������ģʽѡ��
%model = 1;                                                                
if model == 0
    y = Global_controller(tao,x,t,v_tao);
else    
    y = Point_controller(tao,x,t,v_tao);
end
[Fxt,Fyt,M_z] = Force_Moments_calculation(tao,x,y);                        %����ϵͳ��������
global u_trim_now;
if isnan(y(1)) || isnan(y(2)) || isnan(Fxt) || isnan(Fyt) || isnan(M_z) || y(2) <= 0 || abs(y(1)) > 20 || abs(Fxt) > 1000 %�쳣���
     y = u_trim_now;
    [Fxt,Fyt,M_z] = Force_Moments_calculation(tao,x,y); 
end
global controller_output; global x_trim_now; controller_output = [tao;y;x_trim_now;tao;u_trim_now];
 
dx1 = (Fxt*cos(x(3)) - Fyt*sin(x(3)))/2;                                   %״̬���̣�����ϵͳ״̬΢��          
dx2 = (Fxt*sin(x(3)) + Fyt*cos(x(3)))/2;
dx3 = x(4);
dx4 = M_z/0.0514;
dx5 = x(2);
sys = [dx1;dx2;dx3;dx4;dx5];


%{
��������mdlDerivatives
���ܣ����� S_Function�����
���ý����sys ����ϵͳ��� ״̬   ���� ӭ�� ������ ���������� �߶�
example:
 %}
function sys = mdlOutputs(t,x)
global controller_output;
if t == 0
    controller_output = [0;0;0;0;0;0;0;0;0;0;0];
end
sys = [x;controller_output];


%{
��������paraComputation_fixedwing
���ܣ����ݵ�ǰϵͳ״̬�Ϳ������룬����ϵͳ��������
 %}
function [Fxt,Fyt,M_z] = Force_Moments_calculation(tao,x,u)  
tau = tao * pi / 180;                                                      %����ת��
cy_alfa = 0.0247; cy_wz = 2.91;   cy_deltaz = 0.00737;                     %���������������
mz_alfa = -2.04;  mz_wz = -5.865; mz_deltaz = -0.0184;
cy_range = 1.0;%0.7;
mz_range = 1.0;
cy_alfa = cy_range*cy_alfa; cy_wz = cy_range*cy_wz; cy_deltaz = cy_range*cy_deltaz; %�����㶯
mz_alfa = mz_range*mz_alfa; mz_wz = mz_range*mz_wz; mz_deltaz = mz_range*mz_deltaz;

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


%{
��������Global_controller
���ܣ���ƿ��������������
 %}
function y = Global_controller(tao,u,t,v_tao)
x = u(1:5);

p1 = 7.7191e-08; p2 = 6.8331e-06; p3 = -0.00018488; p4 = 0.0012743; p5 = -0.0034819; p6 = 0.02417; p7 = 0.044811;
theta_refer = p1 * t^6 + p2 * t^5 +  p3 * t^4 + p4 * t^3 + p5 * t^2 + p6 * t + p7;
global x_trim_online; x_trim_online = [v_tao;0;theta_refer;0;20];               %ָ����ʼ��ƽ��
global u_trim_online; u_trim_online = [0.3880;0.1400];
global last_t;
if t == 0
    last_t = 0;
end
%if t - last_t > 0.5 || t == 0                                              %��������
[x_trim,u_trim,K_trim,tmin,r0,ce0] = Robust_Stabilize_control(x_trim_online,u_trim_online);
    global K_online; K_online = K_trim;
    global x_trim_now;x_trim_now = x_trim;
    global u_trim_now; u_trim_now = u_trim;
    last_t = t;
%end
    
y = K_online*(x - x_trim_now) + u_trim_now;                                %��ϵͳ�ܿ�������
  

%{
��������system_corridor
���ܣ�������ת��-���� ��������
 %}
function [tao,v_tao] = systemcorridor(t)
t_wing = 0; t_hover = 0;
global t_tran;t_tran = 12.2474;                                            %����ʱ��
global taov0;taov0 = 75; global v0_down;v0_down = 6;                       %�ն�����
v0_up = 23.5; v_tran = v0_up - v0_down; %78������ 


if t < t_wing
    tao = 0;
    v_tao = -(23.5 - v0_up) / t_wing * t + 23.5; 
elseif t > (t_wing + t_tran)
    tao = taov0;
    v_tao = -(v0_down - 3) / t_hover * (t - (t_wing + t_tran)) + v0_down; 
else 
    if t < (t_wing + t_tran/2)                                             %��ת������1
       % tao =  1 * (t - t_wing) * (t - t_wing);
    else
       % tao = -1 * (t  - t_wing - 12.2474) * (t  - t_wing - 12.2474) + 75;
    end
    
    if t == 0
        tao = 0;
        v_tao = v0_up;
    elseif t >= t_tran
        tao = taov0;
        v_tao = v0_down;
    else
        tao = -taov0 * (cos(pi * t / t_tran) + 1) / 2 + taov0;             %��ת��
        
        v_tao = (acos(2 .* tao ./ taov0 - 1)) .* v_tran ./ pi + v0_down;    %�������� S1
        %v_tao =  -8.7e-7.*power(tao,4) + 0.00010316 .* power(tao,3) - 0.005553 .*power(tao,2) - 0.044151 .* tao + 23.43 + 0.03; %��������S Old
        
        
        %v_tao = -v_tran / taov0 .* tao + v0_up;                           %�������� S2
        
        
        %p1 = -2.3486e-09; p2 = 3.9863e-07; p3 = -2.5573e-05; p4 = 0.00076354; p5 = -0.010387; p6 = -7.8404e-05; p7 = 23.499;
        %v_tao = p1*tao^6 + p2*tao^5 + p3*tao^4 + p4*tao^3 + p5*tao^2 + p6*tao + p7;%�������� S3
    end
end

%{
if t < 62.45
    test = 0.01 * t * t;
else
    test = -0.01 * (t - 124.9000) * (t - 124.9000) + 78;
    %{
    if t < 82.4736
        test = -0.01 * (t - 124.9000) * (t - 124.9000) + 78;
    elseif t >= 82.4736 && t < 92.4736
        test = 60;
    else
        test = -0.01 * (t - 10 - 124.9000) * (t - 10 - 124.9000) + 78;
    end
    %}
end
%}