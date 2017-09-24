%���ܣ�����������ѧϵͳģ�ͣ�����tao��V�����е�����ƽ�����Ի������������
%�����ʵ�ֵ�ǰ��������򶨿���
%ע�ǣ���ƽԭ��
%      ���Ի�����ı�����
%      ����Ӧ�õĿ���˼�������
%{
��������Robust_Stabilize_control
���ܣ���ⵥ�������
example:
 %}
function [x_trim,u_trim,K_trim,tmin,r0,ce0,A_trim] = Robust_Stabilize_control(x_trim_set,u_trim_set)
[x_trim,u_trim] = Point_trim(x_trim_set,u_trim_set);                       %������ϵͳ���ж��׶���ƽ
[A_trim,B_trim,~,~] = linmod('trim_realmodel_simulation',x_trim,u_trim);   %����ʵ���ϵͳ������ƽ�����Ի�
%[K_trim,tmin,r0,ce0] = Robust_H_infinite(A_trim,B_trim);                  %³��H���������
[K_trim,tmin,r0,ce0] = Robust_Stabilize(A_trim,B_trim);                    %³���򶨿�����  
end


%{
��������Point_trim
���ܣ���ƽ
example:
 %}
function [x_trim,u_trim] = Point_trim(x_trim_set,u_trim_set)
%requests = [1;2;5];                                                        %��һ�׶���ƽ��Ҫ���ٶȣ��߶ȣ������ǣ�������䵼����
requests = [1;2;3;5];                                                      %�ڶ��׶���ƽ��Ҫ���ٶȣ��������Լ��߶ȣ����������ʺ����룩
y_trim_set = x_trim_set;
global Trimfailed_flag;Trimfailed_flag = 0;

[x_trim,u_trim,~,~] = trim('trim_virtualmodel_simulation',x_trim_set,u_trim_set,y_trim_set,requests,[],requests); %���õ���ƽ�Ƿ�ʧ��
global u_trimlast;global x_trimlast;
if Trimfailed_flag == 1               
    x_trim = x_trimlast;
    u_trim = u_trimlast;
    x_trim(1) = x_trim_set(1);
    Trimfailed_flag = 0;
end

global velocity_up_count;
global t_online; t = t_online;
global flag_velocity;
if t <= 1
    velocity_up_count = 0;
end
if t > 1 && (x_trim(1) - x_trimlast(1)) > 0
    velocity_up_count = velocity_up_count + 1;
    if velocity_up_count == 1
        flag_velocity =  x_trimlast(1);
    end
end

global v_tao_online; 
if abs(x_trim(1) - v_tao_online) > 1 ||  u_trim(1) > 20 || u_trim(1) < -20 || velocity_up_count > 2  %���õ���ƽ�Ƿ��쳣
    fprintf('----------------------------------------------------------------------\n');
    fprintf('----------------------------------------------------------------------\n');
    fprintf('                      Trim abnormal situation                         \n');
    abs(x_trim(1) - x_trimlast(1))
    abs(u_trim(1) - u_trimlast(1))
    x_trim = x_trimlast;
    u_trim = u_trimlast;
    %x_trim(1) = flag_velocity;
    x_trim(1) = x_trim_set(1);
    velocity_up_count = 0;
end
%x_trim = x_trim_set;
x_trimlast = x_trim;
u_trimlast = u_trim;
end


function [K_trim,tmin,r0,ce0] = Robust_H_infinite(A_trim,B_trim)
B1 = B_trim;
C1 = A_trim;
D12 = B_trim;
D11 = B1;
z_define = zeros(5,5); z_define(1,1) = 0; z_define(4,4) = 0; z_define(5,5) = 1;
C1 = z_define * C1; D12 = z_define * D12; D11 = z_define * D11;
%C1 = 0;D12 = 0;D11 = 0;
%D11 = 0;

%LMI��ʼ����LMIϵͳ����
setlmis([]);
X = lmivar(1,[5 1]);
W = lmivar(2,[2 5]);
r = lmivar(1,[1 1]);
%%r = 1;
ce = 7.5;%7.5;
%LMI��������
%LMI_1
lmiterm([1 1 1 X],A_trim,1,'S');
lmiterm([1 1 1 W],B_trim,1);
lmiterm([1 1 1 -W],1,B_trim');
lmiterm([1 2 1 0],B1');
lmiterm([1 2 2 0],-1);

%%{
%ȷ��ϵͳ�������
lmiterm([1 3 1 X],C1,1);
lmiterm([1 3 1 W],D12,1);

lmiterm([1 3 2 0],D11);

%lmiterm([1 3 3 0],-r^2);
lmiterm([1 3 3 r],-1,1);
%ͬʱ����˵��ԭ��ϵͳ�Ĳ�ȷ���Իᵼ�¿��������ɽ�
%}

%LMI_2
lmiterm([-2 1 1 X],1,1);
%LMI_3
lmiterm([3 1 1 r],-1,1);
lmiterm([3 1 1 0],0.1);
%LMI_4
%lmiterm([4 1 1 r],1,1);
%lmiterm([4 1 1 0],-200);
%}
%���LMIϵͳ���ڲ�����
lmisys = getlmis;

%%{
%����feasp��LMIS�Ŀ�����
options = [0;100;1e09;10;1];
[tmin,xfeas] = feasp(lmisys,options);
X0 = dec2mat(lmisys,xfeas,X);
W0 = dec2mat(lmisys,xfeas,W);
r0 = dec2mat(lmisys,xfeas,r);
%r0 = r;
ce0 = ce;
%}

K_trim = W0/X0;%W0 * inv(X0);
end


function [K_trim,tmin,r0,ce0] = Robust_Stabilize(A_trim,B_trim)
load('Matrix3.mat');                                  %����ϵͳģ�͵ĸ���֪��������
alfa_stable = 0.5; 
%Fb = 10 * Fb;

%LMI��ʼ����LMIϵͳ����
setlmis([]);
X = lmivar(1,[5 1]);
W = lmivar(2,[2 5]);
r = 1;
%ce = lmivar(1,[1 1]);
ce = 7.5;%7.5;
%LMI��������
%LMI_1
lmiterm([1 1 1 X],A_trim,1,'S');
lmiterm([1 1 1 W],B_trim,1);
lmiterm([1 1 1 -W],1,B_trim');
lmiterm([1 1 1 X],2 * alfa_stable,1);
%%{
%��ȷ��ϵͳ�������
lmiterm([1 2 1 0],ce * E');
lmiterm([1 2 2 0],-ce);
lmiterm([1 3 1 X],Fa,1);
lmiterm([1 3 1 W],Fb,1);
lmiterm([1 3 3 0],-ce);
%}

%{
%��ȷ��ϵͳ�������
lmiterm([1 2 1 ce],E',1);
lmiterm([1 2 2 ce],-1,1);
lmiterm([1 3 1 X],Fa,1);
lmiterm([1 3 1 W],Fb,1);
lmiterm([1 3 3 ce],-1,1);
%}

%LMI_2
lmiterm([-2 1 1 X],1,1);
%LMI_3
%lmiterm([-3 1 1 ce],1,1);
%}
%���LMIϵͳ���ڲ�����
lmisys = getlmis;

%%{
%����feasp��LMIS�Ŀ�����
options = [0;100;1e09;10;1];
[tmin,xfeas] = feasp(lmisys,options);
X0 = dec2mat(lmisys,xfeas,X);
W0 = dec2mat(lmisys,xfeas,W);
r0 = r;
ce0 = ce;
%}

K_trim = W0/X0;%W0 * inv(X0);
end