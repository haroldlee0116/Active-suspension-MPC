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

%% Discrete the plant 
ts=0.01;          % sampling time
P_dis=c2d(P,ts);  % discrete system P
P_dis.InputName = {'F_{act}';'dot_zr'};
P_dis.OutputName = {'a_s';'(zs - zu)';'zt';};
LTI.A=P_dis.A;
LTI.B=P_dis.B;
LTI.C=P_dis.C;
LTI.D=P_dis.D;

% Definition of system dimension
dim.nx = size(LTI.A,2);     % state dimension
dim.nu = size(LTI.B,2);     % input dimension
dim.ny = size(LTI.D,1);     % output dimension
dim.N = 4;                  % prediction horizon

% Controllability check: rank(A,B)=full row rank
if rank(ctrb(LTI.A,LTI.B))==dim.nx
   disp('Discrete system is controllable')
end
% Observability check: Number of unobservable states
if length(A)-rank(obsv(LTI.A,LTI.C))==0
    disp('Discrete system is observable')
end

%% Definition of quadratic cost function
% stage cost weight
weight.Q = diag([.1,5,.1,5]);                    %weight on state
weight.R = eye(dim.nu)*0.0000001;               %weight on input

% state & input upper and lower bounds
xlb = [-0.088;-0.163;-0.015;-1.965]; %state lower bound
xub = [0.105;0.14;0.029;2.78];       %state upper bound

ulb = [-2500;-inf()];                %input lower bound
uub = [2500; inf()];                 %input upper bound 

% Find LQR.
[K, P] = dlqr(LTI.A, LTI.B, weight.Q, weight.R);
%calculates the optimal gain matrix K such that the state-feedback law
%minimizes the quadratic cost function; P is the infinite horizon solution 
%of the associated discrete-time Riccati equation
K = -K; % Sign convention.
weight.P=P;
weight.QK=weight.Q+K'*weight.R*K;

%% Bump road generation
V=40;              % Speed of vehicle (unit: km/h)
V=V*(1000/3600);   % Speed of vehicle (unit: m/s)
time_length=3;     % Time_length of the nosie (unit: s)
time_sample=ts;    % Sampling time of the noise (unit: s)
bump_L=15;         % The length of the bump(unit: m)
bump_A=0.2;        % The height of the bump(unit: m)
t0=0.8;            % The initial time of the bump(unit:s)
road1 = bumproad(time_length,time_sample,t0,V,bump_L,bump_A);
% The output of the bumproad is dot_zr 
% (the vertical speed of the generated road with a bump (unit:m/s))

%% Compute X_f.
% X_f based on algorithm 1 of exercise set 4 question 2
Xn = struct();
V = struct();
Z = struct();

[Xn.('lqr'), V.('lqr'), Z.('lqr')] = findXn(LTI.A, LTI.B, K, dim.N, xlb, xub, ulb, uub, 'lqr');
terminal = Xn.lqr{1}; % LQR terminal set Xf.

% Post-processing of X_f : 
% Eliminate zero-rows
nonzero_index = [];
for i = 1:size(terminal.A,1)
    
    if all(terminal.A(i,:)==0) == 0
        nonzero_index = [nonzero_index i];
    end
    
end
terminal.A = terminal.A(nonzero_index,:);
terminal.b = terminal.b(nonzero_index);
% Eliminate repeated rows
norep = unique([terminal.A terminal.b], 'rows');

% X_f based on algorithm 1 of exercise set 4 question 2
term.A = norep(:,1:end-1);
term.b = norep(:,end);
% check the initial states in X_f or not
LTI.x0 = [0.05;-0.05;0.015;2];
if ~isempty(find(term.A*LTI.x0-term.b>0))
    disp('The initial state is outside X_f but within X_N')
else 
    disp('The initial state is within X_f')
end
%% Regulation MPC
% Generation of prediction model 
predmod=predmodgen_stability(LTI,dim);            
[H,h]=costgen_stability(predmod,weight,dim);

%simulation horizon
T=300-dim.N; 

%Receding horizon implementation
x=zeros(dim.nx,T+1);
u_rec=zeros(dim.nu,T);
x(:,1)=LTI.x0;

% Input constriant for the Acturator Force: 
b_force = ones(dim.nu*(dim.N),1)*2500;
A_force = [eye(dim.N); -eye(dim.N)];

% State constraint
A_state = [eye(dim.nx);-eye(dim.nx)];
b_state = [0.105;0.14;0.029;2.78;0.088;0.163;0.015;1.9651];

% Terminal constraint: X(N) <= X_f
% term.A*X(N)<=term.b
% X(N)=T_N*x0+S_N*u_N
T_N = predmod.T(end-dim.nx+1:end,:);
S_N = predmod.S(end-dim.nx+1:end,:);
options = sdpsettings('verbose',0,'solver','quadprog');
%% MPC optimization iteration
for k=1:T
    
    x_0=x(:,k);
    
    % Solve the unconstrained optimization problem (with YALMIP)
    u_con = sdpvar(dim.nu*dim.N,1);                    % define optimization variable
                                               
                                        
    Constraint=[u_con(2:2:end)==road1(k:(dim.N)-1+k)'; % Inputdot(z_r) following the road profile
                A_force*u_con(1:2:end)<=b_force;       % Acturator force boundary
                A_state*x_0<=b_state;                  % State constraints
                term.A*(S_N*u_con+T_N*x_0)<=term.b;]; % Terminal set

    Objective = 0.5*u_con'*H*u_con+(h*x_0)'*u_con;     % define cost function
    optimize(Constraint,Objective,options);            % solve the problem
    u_con=value(u_con);                                % assign the solution
    
    % Select the first input only
    u_rec(:,k)=u_con(1:dim.nu);

    % Compute the state/output evolution
    x(:,k+1)=LTI.A*x_0 + LTI.B*u_rec(:,k);
    clear u_con
    
end
% save the state trajectories and input sequence for plots
x_outsidexf3=x;
u_rec_outsidexf3=u_rec;
save x_outsidexf3
save u_rec_outsidexf3
%% Plots state trajecotries and input sequence for three cases
 
load('u_rec_insidexf.mat')  %initial state inside X_f
load('u_rec_outsidexf2.mat')%initial state outside X_f but inside X_N
load('u_rec_outsidexf3.mat')%initial state inide X_N + bump road disturbance
load('x_insidexf.mat')
load('x_outsidexf2.mat')
load('x_outsidexf3.mat')
 
        figure(1)
        subplot(2,3,1)
        hold on
        plot(0:T, x_insidexf(1,:));
        plot(0:T, x_outsidexf2(1,:));
        plot(0:T, x_outsidexf3(1,:));
        title('State: Suspension Stroke')
        ylabel('z_s-z_u (m)')
        xlabel('simulation steps (0.01s/step)')
        legend('In X_f without bump road','In X_N without bump road','In X_N with bump road')
        hold off
        
        subplot(2,3,2)
        hold on
        plot(0:T, x_insidexf(2,:));
        plot(0:T, x_outsidexf2(2,:));
        plot(0:T, x_outsidexf3(2,:));
        ylabel('dot(z_s) (m/s)')
        xlabel('simulation steps (0.01s/step)')
        title('State: Sprung Mass Velocity')
        hold off
        
        subplot(2,3,3)
        hold on
        plot(0:T, x_insidexf(3,:));
        plot(0:T, x_outsidexf2(3,:));
        plot(0:T, x_outsidexf3(3,:));
        title('State: Tyre Deflection')
        ylabel('zt (m)')
        xlabel('simulation steps (0.01s/step)')
        hold off
        
        subplot(2,3,4) 
        hold on
        plot(0:T, x_insidexf(4,:));
        plot(0:T, x_outsidexf2(4,:));
        plot(0:T, x_outsidexf3(4,:));
        title('State: Unsprung Mass Velocity')
        ylabel('dot(z_u) (m/s)')
        xlabel('simulation steps (0.01s/step)')
        hold off
        
        subplot(2,3,5)
        hold on
        plot(1:T, u_rec_insidexf(1,:));
        plot(1:T, u_rec_outsidexf2(1,:));
        plot(1:T, u_rec_outsidexf3(1,:));
        title('Input: Actuator Force (N)')
        ylabel('Actuator Force (N)')
        xlabel('simulation steps (0.01s/step)')
        hold off
        
        subplot(2,3,6)
        hold on
        plot(1:T, u_rec_insidexf(2,:));
        plot(1:T, u_rec_outsidexf2(2,:));
        plot(1:T, u_rec_outsidexf3(2,:));
        title("Bump Road Preview Disturbance in Regulation MPC")
        ylabel('z_s (m/s)')
        xlabel('simulation steps (0.01s/step)')
        hold off
%% CLF decrease
Vf=zeros(T-1,1);
Vfplus=zeros(T-1,1);
LHS=zeros(T-1,1);
l=zeros(T-1,1);
for i=1:T-1
    Vf(i) =0.5*x(:,i)'*weight.QK*x(:,i)+0.5*x(:,i)'*weight.P*x(:,i);
    Vfplus(i) = 0.5*x(:,i+1)'*weight.QK*x(:,i+1)+0.5*x(:,i+1)'*weight.P*x(:,i+1);
    LHS(i) = Vfplus(i)-Vf(i);
    l(i) = 0.5*x(:,i)'*weight.QK*x(:,i);
end

for i=1:T-1
    Vf(i) =0.5*x(:,i)'*weight.P*x(:,i);
    Vfplus(i) = 0.5*x(:,i+1)'*weight.P*x(:,i+1);
    LHS(i) = Vfplus(i)-Vf(i);
    l(i) = 0.5*x(:,i)'*weight.QK*x(:,i);
end
figure(3)
plot(LHS(1:50),'--b','LineWidth',1)
hold on
plot(-l(1:50),'r:','LineWidth',1)
legend('V_f(x+1)-V_f(x)','-l(x,u)')
xlabel('simulation horizon')
grid on
hold off
%% X_N region of attraction
% Estimate of the region of attraction with the algorithm synthesized in
% exercise "Computing the projection set"

% Prepare H_u
H_u=zeros(dim.N,dim.N*dim.nu);
column=1;
for i=1:dim.N
    H_u(i,column)=1;
    column=column+2;
end
H_u=[H_u;-H_u];

% Prepare Psi_state (manuelly)
Psi_state_upper=[b_state(1:4);b_state(1:4);b_state(1:4);b_state(1:4)];
% the state should within the upper bound, size=(dim.nx*N,1)
Psi_state_lower=[b_state(5:8);b_state(5:8);b_state(5:8);b_state(5:8)];
% the state should within the lower bound, size=(dim.nx*N,1)
Psi_state=[Psi_state_upper(); Psi_state_lower];

% Polytope Z : Z=:={(x,u)|G*x_0+H*u_N+Psi<=0}
% [ T|      [ S |        [state upper bound|
% |-T| *x + |-S |*u_N <= |state lower bound| 
% | 0]      |H_u]        |     b_force     ]
G=[predmod.T(1:end-dim.nx,:);-predmod.T(1:end-dim.nx,:);zeros(size(H_u,1),dim.nx)];
H=[predmod.S(1:end-dim.nx,:);-predmod.S(1:end-dim.nx,:);H_u];
Psi=-[Psi_state;b_force];


G_i = [G H(:,1:end-1)];
H_i = H(:,end);
Psi_i = [Psi];

for i = dim.N-1:-1:0

    [P_i, gamma_i] = single_input(G_i,H_i,Psi_i);
    
    G_i = P_i(:,1:end-1);
    H_i = P_i(:,end);
    Psi_i = gamma_i;
    
end

P = P_i;
gamma = gamma_i;


% Post-processing
% Eliminate zero-rows
nonzero_index = [];
for i = 1:size(P,1)
    
    if all(P(i,:)==0) == 0
        nonzero_index = [nonzero_index i];
    end
    
end

P = P(nonzero_index,:);
gamma = gamma(nonzero_index);

% Eliminate repeated rows
norep = unique([P gamma], 'rows');

P = norep(:,1:end-1);
gamma = norep(:,end);

% Plot region
x_1 = -0.15:.01:0.15;
x_2 = -0.19:.01:0.17;

[X_1, X_2] = meshgrid(x_1, x_2);
clear map;

for i = 1:size(P,1)
    map(:,:,i) = (P(i,2).*X_2 <= -P(i,1).*X_1 - gamma(i));    
end

figure(4)
contourf(X_1, X_2, min(double(map), [], 3), [1 1],':');
xlabel('$z_s-z_u $ (state constriants:[-0.08 0.105])','Interpreter','latex'), ylabel('$\dot{zs} $ (state constriants: [-0.163 0.14])','Interpreter','latex'), grid on;
xlim([-0.15 0.15]), ylim([-0.19 0.17]);
title('$x_N$ for different state x1 and x2 initial conditions (Initial state is red circle)','Interpreter','latex')
hold on
plot([1 1]*-0.08, ylim, '--k')               
plot([1 1]*0.105, ylim, '--k')               
plot(xlim,[1 1]*-0.163, '--k')               
plot(xlim,[1 1]*0.14, '--k')
hold on
scatter(LTI.x0(1,:),LTI.x0(2,:),'red')
hold off
