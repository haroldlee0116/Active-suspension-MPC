close all; clear all; clc

%% Plant model
m_s = 395.3;      % kg
m_u = 48.3;       % kg
k_s = 30.01e+3;    % N/m
c_s = 1450;       % Ns/m
k_t = 3.4e5;      % N/m

% Inputs: u = [Fact dot_zr]                   
% States: x = [zs-zu; dot_zs; zt; dot_zu]    
% Output: y= [a_s; susTravel; zt]
% y(1) = a_s = dotdot_zs = [-k_s/m_s -c_s/m_s k_s/m_s c_s/m_s]x + [1/m_s 0]u
% y(2) = susTravel = zs - zu = [1 0 -1 0]x +[0 0]u
% y(3) = zt = zu-zr = [0 0 1 0]x + [0 -1]

% Equation: dot_x = Ax + Bu
% Measurements: y = Cx + Du 

A = [0 1 0 -1;-k_s/m_s -c_s/m_s 0 c_s/m_s;...
    0 0 0 1;k_s/m_u c_s/m_u -k_t/m_u -c_s/m_u];
B = [0 0;-1/m_s 0;0 -1;-1/m_u 0];    
C = [-k_s/m_s -c_s/m_s 0 c_s/m_s;...
    1 0 0 0;...
    0 0 1 0];                             
D = [1/m_s 0;0 0;0 0];

P = ss(A,B,C,D);
P.InputName = {'F_{act}';'dot_zr'};
P.OutputName = {'a_s';'(zs - zu)';'zt';};

% Discrete the plant 
ts=0.01;          % sampling time
P_dis=c2d(P,ts);  % discrete system P
P_dis.InputName = {'F_{act}';'dot_zr'};
P_dis.OutputName = {'a_s';'(zs - zu)';'zt';};
LTI.A = P_dis.A;
LTI.B = P_dis.B;
LTI.C = P_dis.C;
LTI.D = P_dis.D;
x0 = [0; 0; 0; 0]; %Initial state
%% Bump road generation
V=40;              % Speed of vehicle (unit: km/h)
V=V*(1000/3600);   % Speed of vehicle (unit: m/s)
time_length=3;     % Time_length of the nosie (unit: s)
time_sample=ts;    % Sampling time of the noise (unit: s)
bump_L=15;         % The length of the bump(unit: m)
bump_A=0.4;        % The height of the bump(unit: m)
t0=0.5;             % The initial time of the bump(unit:s)
bumproad1= bumproad(time_length,time_sample,t0,V,bump_L,bump_A);
%The output of the bumproad is dot_zr 
%(the vertical speed of the generated road with a bump (unit:m/s))
%% LTI system definition
% x+ = A*x + B*u + Bd*d
% y  = C*x + D*u + Cd*d
LTI.Cd=[ 0.001 0; 0 0.01; 0.00001 0];
LTI.Bd=[ 0.01 0; 0 0.01; 0 0.001; 0 0.01]; 
LTI.x0=x0;
LTI.d=[0.1; 0.1];
LTI.yref=[0; 0; 0];

%Definition of system dimension
dim.nx = size(LTI.A,2);     % state dimension
dim.nu = size(LTI.B,2);     % input dimension
dim.ny = size(LTI.D,1);
dim.nd = 2;          %disturbance dimension
dim.N = 5;           %horizon

%% Definition of quadratic cost function
weight.Q = diag([1 1 1 1])*1;                        %weight on output
weight.R = eye(dim.nu)*0.0001;                       %weight on input
weight.P = dare(LTI.A,LTI.B,weight.Q,weight.R);      %terminal cost
T=300;     %simulation horizon

% Check if the problem is well posed
rank_check = rank([eye(dim.nx)-LTI.A -LTI.Bd; LTI.C LTI.Cd]);
rank_shouldbe = dim.nx+dim.nd;

if rank_check == rank_shouldbe
    if dim.nd <= dim.ny
        disp("The augemented system is observable, and nd<=p")
    end
end

%% Extended system computation

LTIe.A=[LTI.A LTI.Bd; zeros(dim.nd,dim.nx) eye(dim.nd)];
LTIe.B=[LTI.B; zeros(dim.nd,dim.nu)];
LTIe.C=[LTI.C LTI.Cd];
LTIe.D=LTI.D;   
LTIe.x0=[LTI.x0; LTI.d];
LTIe.yref=LTI.yref;

%Definition of system dimension
dime.nx=dim.nx+dim.nd;     %state dimension
dime.nu=dim.nu;            %input dimension
dime.ny=dim.ny;            %output dimension
dime.N=5;                  %horizon


%Definition of quadratic cost function for extened system
weighte.Q=blkdiag(weight.Q,zeros(dim.nd));            %weight on output
weighte.R=weight.R;                                   %weight on input
weighte.P=blkdiag(weight.P,zeros(dim.nd));            %terminal cost
%% Offset-free MPC from output ：Observe output, not observe x & d
 
predmode=predmodgen_outputMPC(LTIe,dime); 
[He,he]=costgen_outputMPC(predmode,weighte,dime); 

% Receding horizon implementation
xe=zeros(dime.nx,T+1);
y=zeros(dime.ny,T+1);
u_rec=zeros(dime.nu,T);
xehat=zeros(dime.nx,T+1);

xe(:,1)=LTIe.x0;
xehat(:,1)=[0.1; 0; 0; 0.5; 0; 0];
y(:,1)=LTIe.C*LTIe.x0;

% Simulation horizon
steps=T-dime.N+1;  
% Input constriant for the Acturator Force: A*u_uncon(1:2:end)<=b
b_force = ones(dime.nu*(dime.N),1)*2500;
A_force = [eye(dime.N); -eye(dime.N)];
% Observer gain
L=place(LTIe.A',LTIe.C',[0.5; 0.4; 0.45;0.6;0.65;0.8])';
xr_record=zeros(dim.nx,steps);
ur_record=zeros(dim.nu,steps);
options = sdpsettings('verbose',0,'solver','quadprog');
for k=1:steps
    
    xe_0=xe(:,k);  
    dhat=xehat(end-dim.nd+1:end,k);
    
    %Compute optimal ss (online, at every iteration)
    eqconstraints=eqconstraintsgen(LTI,dim,dhat);
    [xr,ur]=optimalss(LTI,dim,weight,[],eqconstraints); 
    xre=[xr;dhat];

    uostar = sdpvar(dime.nu*dime.N,1);                      %define optimization variable
    Constraint=[uostar(2:2:end)==bumproad1(k:(dim.N)-1+k)'; %follow the road profile
                A_force*uostar(1:2:end)<=b_force;];         %input constraints
                                                                   
    Objective = 0.5*uostar'*He*uostar+(he*[xe_0; xre; ur])'*uostar;    %define cost function
    optimize(Constraint,Objective,options);                            %solve the problem
    uostar=value(uostar);      

    % Select the first input only
    u_rec(:,k)=uostar(1:dim.nu);

    % Compute the state/output evolution
    xe(:,k+1)=LTIe.A*xe_0 + LTIe.B*u_rec(:,k);
    y(:,k+1)=LTIe.C*xe(:,k+1) + LTIe.D*u_rec(:,k+1);
    clear uostar 
        
    % Update extended-state estimation
    xehat(:,k+1)=LTIe.A*xehat(:,k)+LTIe.B*u_rec(:,k)+L*(y(:,k)-LTIe.C*xehat(:,k));
    % Record the state reference
    xr_record(:,k)=xr;
    ur_record(:,k)=ur;
end
 

%% Estimation and real state/input trajectories plots
%       Real & Predicted input and road disturbance plot
        figure(1)
        subplot(1,2,1)
        hold on
        plot(1:steps, u_rec(1,1:steps));
        plot(1:steps, ur_record(1,1:steps));
        legend('Real input value','Reference from OTS')
        title('Input: Actuator Force (N)')
        ylabel('F_act (N)')
        xlabel('simulation steps (0.01s/step)')
        hold off
        
        subplot(1,2,2)
        hold on
        plot(1:steps, u_rec(2,1:steps),'-o');
        plot(1:steps, ur_record(2,1:steps));
        legend('Bump road profile ','Reference from OTS')
        title('Input: Bump Road Profile (m/s)')
        ylabel('z_r (m/s)')
        xlabel('simulation steps (0.01s/step)')
        hold off

%       Real & Estimated state plot
        figure(2)
        subplot(2,2,1)
        hold on
        plot(1:steps, xe(3,1:steps));
        plot(1:steps, xehat(3,1:steps));
        plot(1:steps, xr_record(3,1:steps));
        legend('Real state value','Estimated state value','Reference')
        title('State: Tyre Deflection')
        ylabel('z_t (m)')
        xlabel('simulation steps (0.01s/step)')
        hold off
        
        subplot(2,2,2)
        hold on
        plot(1:steps, xe(4,1:steps));
        plot(1:steps, xehat(4,1:steps));
        plot(1:steps, xr_record(4,1:steps));
        legend('Real state value','Estimated state value','Reference')
        title('State: Unsprung Mass Velocity')
        ylabel('dot(z_u) (m/s)')
        xlabel('simulation steps (0.01s/step)')
        hold off
        
        subplot(2,2,3)
        hold on
        plot(1:steps, xe(5,1:steps));
        plot(1:steps, xehat(5,1:steps));
        legend('Real state value d1','Estimated state value d1')
        title('State: Disturbance d1')
        ylabel('d1')
        xlabel('simulation steps (0.01s/step)')
        hold off
        
        subplot(2,2,4)
        hold on
        plot(1:steps, xe(6,1:steps));
        plot(1:steps, xehat(6,1:steps));
        legend('Real state value d2','Estimated state value d2')
        title('State: Disturbance d2')
        ylabel('d2')
        xlabel('simulation steps (0.01s/step)')
        hold off
        
%% Discrete system input
% Inputs: [Fact dot_zr]
time_length1=3;
t = time_sample:ts:time_length1;       % sample data length
% Fact
F=ones(length(bumproad1),1)*5;         % size: sample data length * 1
F=u_rec(1,:);                          % size: sample data length * 1
F2=[t'  F']; 
% First column: sample time ; Second column: Fact at sample time

% dot_zr
road2=[t' bumproad1(1:T)'];
% First column: sample time ; Second column: dot_zr at sample time

% Simulink output plots
sim('carModelMPC');

% Outputs: [a_s (zs - zu) zt]
t1 = 0:ts:time_length1;        %Starting from 0 => sample data length+1
output_as=output.data(:,1);    %a_s   (m/s^2)
output_zszu=output.data(:,2);  %zs-zu (m)
output_zt=output.data(:,3);    %zt    (m)
% Plot outputs
figure(3)
plot(t1,output_as); 
legend('a_s');
title('Sprung mass acceleration (m/s^2)');
xlabel('Time (s)');
ylabel('Sprung mass acceleration (m/s^2)');
figure(4)
plot(t1,output_zszu);
legend('zs - zu');
title('Suspension stroke (m)');
xlabel('Time (s)');
ylabel('Suspension stroke (m)');
figure(5)
plot(t1,output_zt);
legend('zt');
title('Tire deflection (m)');
xlabel('Time (s)');
ylabel('Tire deflection (m)');

        
      