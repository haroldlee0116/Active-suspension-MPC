close all; clear all; clc
% In this file, the velocity update algorithm is adopted to find the
% maximum velocity that makes the optimization problem feasible while
% meeting all the constraints

% Note: The plots in the report are not from this code. Because of the paper
% size requirment, we feel very sorry to not provide the results of this
% algorithm in the report. But we feel very happy to work out this
% algorithm that is useful in the practice. 
%% Plant model
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
Ad = P_dis.A;
Bd = P_dis.B;
Cd = P_dis.C;
Dd = P_dis.D;
%% Road generation
% Variables
road_class=3;   % Values For ISO Road C Roughness Classification, from 1 to 8
V=30;           % Speed of vehicle (unit: km/h)
V=V*(1000/3600);% Speed of vehicle (unit: m/s)
time_length=3;  % Time_length of the nosie (unit: s)
time_sample=ts; % Sampling time of the noise (unit: s)

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
x0 = [0 ;0 ;0;0];

% Definition of system dimension
dim.nx = size(LTI.A,2);     % state dimension
dim.nu = size(LTI.B,2);     % input dimension
dim.ny = size(LTI.D,1);
dim.N = 3;                  % prediction horizon

%Definition of quadratic cost function
Q = [0.01 0 0; 0 .00001 0; 0 0 .000001];   % weight on output
R = eye(dim.nu)*0.0000000000001;           % weight on input
 
% Generation of prediction model
[P_state,S_state]=predmodgen_state(LTI,dim);    

% Set some options for YALMIP and solver
options = sdpsettings('verbose',0,'solver','quadprog');

% Receding horizon implementation for the unconstrained control problem
T=296;
x_0 = x0;
x(:,1) = x0;

% Input constriant for the Acturator Force: A*u_uncon(1:2:end)<=b
b_force = ones(dim.nu*(dim.N),1)*2500;
A_force = [eye(dim.N); -eye(dim.N)];

% Output constraint
A_output = [1 0 0;-1 0 0; 0 1 0; 0 -1 0; 0 0 1];
b_output = [1;1; 0.1; 0.1 ; 0.01];

%% Velocity update optimization
disp('Original Velocity is (m/s)'),disp(V)
%The iteration number for velocity update algorithm
alpha = [1 1 1 1 1 1 1 1];   
for i=1:length(alpha) 
    for k=1:T

        % Write the cost function in quadratic form
        [H,h,const]=costgen(P_state,S_state,alpha(i)*Q,R,dim,x_0); 
        
        % Solve the constrained optimization problem (with YALMIP)
        u_uncon = sdpvar(dim.nu*dim.N,1); % define optimization variable

        Constraint=[u_uncon(2:2:end)==road1(k:(dim.N)-1+k)';
                     A_force*u_uncon(1:2:end)<=b_force;
                     A_output*(LTI.C*x_0+LTI.D*u_uncon(1:2))<=b_output;];
                    % (1) first equality constraint:
                    % to control the second input equal to the road profile 
                    % (2) second inequality constraint:
                    % to control the first input within the boundary
                    % (3) third inequality constraint:
                    % to control the output within the comfortable and safe
                    % boundary
        Objective = 0.5*u_uncon'*H*u_uncon+h'*u_uncon;  %define cost function
    
        %solve the problem
        diagnostics = optimize(Constraint,Objective,options);
        if diagnostics.problem == 0
             if k==T
                 V=V+3;
                 disp('Solver thinks it is feasible, velocity has been used for the plots. Velocity increase to (m/s)')
                 disp(V)
             end
        elseif diagnostics.problem == 1
                V=V-2;
                road1=roadGenerator2(road_class,V,time_length,time_sample);
                disp('Solver thinks it is infeasible, Velocity decease to (m/s)')
                disp(V)
                clear u_uncon u_rec x
                x_0 = x0;
                break
        end
        u_uncon=value(u_uncon);                 

        % Select the first input only
        u_rec(:,k) = u_uncon(1:dim.nu);
        
        % Compute the state/output evolution
        x(:,k+1)=LTI.A*x_0 + LTI.B*u_rec(:,k);
        
        % Update initial state for the next iteration
        x_0=x(:,k+1);

        clear u_uncon；
    end
    if V<=0 
        disp('V has reached zero, but the problem is still not feasible')
        break
    end
    
    if diagnostics.problem == 0
        figure(1)
        subplot(2,2,1)
        plot(0:T, x(1,:));
        legend('z_s-z_u (m)')
        title('the state trajectories: z_s-z_u')        
        subplot(2,2,2)
        plot(0:T, x(2,:));
        legend('dot(z_s) (m/s)')        
        title('the state trajectories: dot(z_s)')      
        subplot(2,2,3)
        plot(0:T, x(3,:));
        legend('zt (m)')
        title('the state trajectories: zt')
        subplot(2,2,4)
        plot(0:T, x(4,:));
        legend('dot(z_u) (m/s)')
        title('the state trajectories: dot(z_u)')
        
        figure(2)
        subplot(1,2,1)
        plot(0:T-1, u_rec(1,:), '-o');
        legend('Fact')
        title('Acuatorc Force')
        subplot(1,2,2)
        plot(0:T-1, u_rec(2,:), '-o');
        legend('z_r')
        title('Road profile')
        x_0 = x0;
        x(:,1) = x0;
        if i ~= length(alpha)
            u_feasible=u_rec;
            road_feasible=road1;
            V_feasible=V-3;
            road1=roadGenerator2(road_class,V,time_length,time_sample);
            clear u_rec 
        else
             u_feasible=u_rec;
             road_feasible=road1; 
            disp('The final feasible velocity is (m/s)')
            disp(V-3)
        end
    end
    if diagnostics.problem == 1
        if V<=0 
            disp('V has reached zero, but the problem is still not feasible')
            break
        elseif i==length(alpha)
            disp('Run out the chance! Lower down the original velocity and rerun the code')
            disp('The final velocity is (m/s)')
            disp(V_feasible);
            break
        end
    end

    
 end

%% Discrete system input
% Inputs: [Acuator force Road progile]
time_length1=2.96;
t = time_sample:ts:time_length1;  % sample data length

% Acuator force
F=u_feasible(1,:);                
F2=[t'  F']; 
% First column: sample time ; Second column: Fact at sample time

% Road progile
road2=[t' road_feasible(1:T)'];
% First column: sample time ; Second column: dot_zr at sample time
%% Simulink output plots
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

