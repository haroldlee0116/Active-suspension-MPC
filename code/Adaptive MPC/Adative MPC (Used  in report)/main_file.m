close all; clear all; clc
% If there is an error showing up in line 110, please rerun the code.
% (It is the feasibility error, it may take several trials)

% We solve the problem of potential unfeasibility in 
% "Velocity update Adaptive MPC", but due to the report length, we did not
% show the results of Velocity update MPC. But it is a very successful MPC
%% Plant model
for V=30:30:60
m_s = 395.3;      % kg
m_u = 48.3;       % kg
k_s = 30.01e+3;   % N/m
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

Ad=P_dis.A;
Bd=P_dis.B;
Cd=P_dis.C;
Dd=P_dis.D;
%% Road generation
% Variables
road_class=3;      % Values For ISO Road C Roughness Classification, from 1 to 8
V1=round(V);       % Speed of vehicle (unit: km/h) %% try first V=30, V=60 
V=V*(1000/3600);   % Speed of vehicle (unit: m/s)
time_length=3;     % Time_length of the nosie (unit: s)
time_sample=ts;    % Sampling time of the noise (unit: s)

% Apply roadGenerator2 function to simulink the road 
road1=roadGenerator2(road_class,V,time_length,time_sample);
% The output of the roadGenerator2 is dot_zr 
% (the vertical speed of the generated whit Gaussian noise road (unit:m/s))
%% Time-varying MPC with constriants
% Definition of the LTI system
LTI.A = Ad;
LTI.B = Bd;
LTI.C = Cd;
LTI.D = Dd;

x0 = [0;0;0;0];

% Definition of system dimension
dim.nx = size(LTI.A,2);     % state dimension
dim.nu = size(LTI.B,2);     % input dimension
dim.ny = size(LTI.D,1);
dim.N = 3;                  % prediction horizon

%Definition of quadratic cost function
Q= [0.0001 0 0;0 0 0;0 0 10]; %weight on output
R = eye(dim.nu);        %weight on input
 
% Generation of prediction model
[P_state,S_state]=predmodgen_state(LTI,dim);    

% Set some options for YALMIP and solver
options = sdpsettings('verbose',0,'solver','quadprog');

% Receding horizon implementation for the unconstrained control problem
T=296;
x_0 = x0;
x(:,1) = x0;

% Input constraint for the Acturator Force: A*u_uncon(1:2:end)<=b
b_force = ones(dim.nu*(dim.N),1)*2500;
A_force = [eye(dim.N); -eye(dim.N)];


% Output constraint
A_output = [1 0 0;-1 0 0; 0 1 0; 0 -1 0; 0 0 1;0 0 -1];

if V1 == 30
 b_output = [0.65;0.65;0.09;0.08;0.0128;0.0128]; %V =30 km/hr
end

if V1 == 60
 Q= [0.0001 0 0;0 0 0;0 0 1];              % V=60 km/hr   
 b_output = [1.35;1.35;0.09;0.08;0.0128;0.0128]; %V =60 km/hr
end

 for k=1:T

        % Write the cost function in quadratic form
        [H,h,const]=costgen(P_state,S_state,Q,R,dim,x_0); 
        % Solve the constrained optimization problem (with YALMIP)
        u_uncon = sdpvar(dim.nu*dim.N,1);        % define optimization variable

        Constraint=[u_uncon(2:2:end)==road1(k:(dim.N)-1+k)';
                     A_force*u_uncon(1:2:end)<=b_force;
                     A_output*(LTI.C*x_0+LTI.D*u_uncon(1:2))<=b_output];
                    % (1) first equality constraint:
                    % to control the second input equal to the road profile 
                    % (2) second inequality constraint:
                    % to control the first input within the boundary
                    % (3) third inequality constraint:
                    % to control the output within the comfortable and safe
                    % boundary
        Objective = 0.5*u_uncon'*H*u_uncon+h'*u_uncon;  %define cost function

        optimize(Constraint,Objective,options);

        u_uncon=value(u_uncon);                  %assign the solution to uopt
        % Select the first input only
        u_rec(:,k) = u_uncon(1:dim.nu);
        
        % Compute the state/output evolution
        x(:,k+1)=LTI.A*x_0 + LTI.B*u_rec(:,k);
        
        % Update initial state for the next iteration
        x_0=x(:,k+1);

        clear u_uncon
 end
        figure(1)
        subplot(2,2,1)
        plot(0:T, x(1,:));
        ylabel('z_s-z_u (m)')
        title('Suspension stroke')
        if V1 == 60
        legend('V = 30 km/hr','V = 60 km/hr');
        end
        subplot(2,2,2)
        plot(0:T, x(2,:));
        hold on;
        ylabel('dot(z_s) (m/s)')
        title('Sprung mass velocity')
        if V1 == 60
        legend('V = 30 km/hr','V = 60 km/hr')
        end
        subplot(2,2,3)
        plot(0:T, x(3,:));
        hold on;
        ylabel('zt (m)')
        title('Tyre Deflection')
        if V1 == 60
        legend('V = 30 km/hr','V = 60 km/hr')
        end
        subplot(2,2,4)
        plot(0:T, x(4,:));
        hold on;
        ylabel('dot(z_u) (m/s)')
        xlabel('t (s)')
        if V1 == 60
        legend('V = 30 km/hr','V = 60 km/hr')
        end
        title('Unsprung mass elevation');
        
        figure(2)
        plot(0:T-1, u_rec(1,:), '-o');
        hold on;
        ylabel('f_A (N)');
        if V1 == 60
        legend('V = 30 km/hr','V = 60 km/hr')
        end
        title('Actuator force');

        
        figure(4)
        plot(0:T-1, u_rec(2,:), '-o');
        hold on;
        if V1 == 60
        legend('V = 30 km/hr','V = 60 km/hr')
        end
        title('Road input')
        x_0 = x0;
        x(:,1) = x0;
        
%% Discrete system input
% Inputs: [Acuator force, road profile]
time_length1=2.96;
t = time_sample:ts:time_length1;  % sample data length

F=u_rec(1,:);            
F2=[t'  F']; 
% First column: sample time ; Second column: Fact at sample time

% dot_zr
road2=[t' road1(1:T)'];
% First column: sample time ; Second column: dot_zr at sample time
%% Simulink output plots
sim('carModelMPC');

% Outputs: [a_s (zs - zu) zt]
time_length2=2.96;
t1 = 0:ts:time_length2;        %Starting from 0 => sample data length+1
output_as=output.data(:,1);    %a_s   (m/s^2)
output_zszu=output.data(:,2);  %zs-zu (m)
output_zt=output.data(:,3);    %zt    (m)
if V1 == 30
   RMSV30_1 = rms(output_as);
%    save('RMS301');
end
if V1 == 60
   RMSV60_1 = rms(output_as); 
%    save('RMS601');
end

% Plot outputs
figure(3)
subplot(3,1,1);
plot(t1,output_as);
ylabel('a_s (m/s^2)');
hold on;
    if V1 == 60
       legend('V = 30 km/hr','V = 60 km/hr')
    end
title('Sprung mass acceleration (m/s^2)');

subplot(3,1,2);
plot(t1,output_zszu);
hold on;
ylabel('z_s - z_u')
    if V1 == 60
        legend('V = 30 km/hr','V = 60 km/hr')
    end
title('Suspension stroke (m)');
subplot(3,1,3);
plot(t1,output_zt);
hold on;
    if V1 == 60
        legend('V = 30 km/hr','V = 60 km/hr')
    end
title('Tire deflection (m)');
ylabel('z_t (m)');

%% Only passive suspension (No MPC)
F2=[t' 0*F']; 

sim('carModelMPC');

% Outputs: [a_s (zs - zu) zt]
time_length2=2.96;
t1 = 0:ts:time_length2;        %Starting from 0 => sample data length+1
output_as=output.data(:,1);    %a_s   (m/s^2)
output_zszu=output.data(:,2);  %zs-zu (m)
output_zt=output.data(:,3);    %zt    (m)
if V1 == 30
   RMSV30_2 = rms(output_as); 
   save('RMS30');
end
if V1 == 60
   RMSV60_2 = rms(output_as);
   save('RMS60');
end
% Plot outputs
figure(5)
subplot(3,1,1);
plot(t1,output_as);
ylabel('a_s (m/s^2)');
hold on;
    if V1 == 60
       legend('V = 30 km/hr','V = 60 km/hr')
    end
title('Sprung mass acceleration (m/s^2)');
subplot(3,1,2);
plot(t1,output_zszu);
hold on;
ylabel('z_s - z_u')
    if V1 == 60
        legend('V = 30 km/hr','V = 60 km/hr')
    end
title('Suspension stroke (m)');
subplot(3,1,3);
plot(t1,output_zt);
hold on;
    if V1 == 60
        legend('V = 30 km/hr','V = 60 km/hr')
    end
title('Tire deflection (m)');
ylabel('z_t (m)');

if V1 == 30
    impV30 = abs((RMSV30_2 - RMSV30_1)/RMSV30_2)*100;
    impV30 = impV30*ones(length(output_as));
end

if V1 == 60
    impV60 = abs((RMSV60_2 - RMSV60_1)/RMSV60_2)*100;
    impV60 = impV60*ones(length(output_as));
end
clear all
end 

%% Improvement MPC vs No MPC
load('RMS30');
load('RMS60');
impV30 = abs((RMSV30_2 - RMSV30_1)/RMSV30_2);
disp('Compared to negatuve suspension system, the comfort in active suspension system for V=30 km/h increase'),disp(impV30)
impV60 = abs((RMSV60_2 - RMSV60_1)/RMSV60_2);
disp('Compared to negatuve suspension system, the comfort in active suspension system for V=60 km/h increase'),disp(impV60)