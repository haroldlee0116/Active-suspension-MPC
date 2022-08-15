close all; clear all; clc
%% Load recorded data from regulation MPC "main_file.m"
load('u_rec_000.mat')
load('u_rec_1.mat')
load('u_rec2.mat')
load('u_rec3.mat')
load('u_rec4.mat')
load('u_rec5.mat')
load('u_rec7.mat')
load('u_rec6.mat')
load('u_rec8.mat')
load('u_rec9.mat')
load('u_rec10.mat')
load('u_rec12.mat')
load('u_rec11.mat')
load('x_000.mat')
load('x_1.mat')
load('x2.mat')
load('x3.mat')
load('x4.mat')
load('x5.mat')
load('x6.mat')
load('x7.mat')
load('x9.mat')
load('x8.mat')
load('x10.mat')
load('x11.mat')
load('x12.mat')
%% Q
        figure(1)
        subplot(2,2,1)
        hold on
        plot(x_000(1,:));
        plot(x_1(1,:));
        plot(x2(1,:));
        title('State: Suspension Stroke')
        ylabel('z_s-z_u (m)')
        xlabel('simulation steps (0.01s/step)')
        hold off
        
        subplot(2,2,2)
        hold on
        plot(x_000(2,:));
        plot(x_1(2,:));
        plot(x2(2,:));
        title('State: Sprung Mass Velocity')
        ylabel('dot(z_s) (m/s)')
        xlabel('simulation steps (0.01s/step)')
        hold off
        
        subplot(2,2,3)
        hold on
        plot(x_000(3,:));
        plot(x_1(3,:));
        plot(x2(3,:));
        title('State: Tyre Deflection')
        ylabel('zt (m)')
        xlabel('simulation steps (0.01s/step)')
        hold off
        
        subplot(2,2,4)
        hold on
        plot(x_000(4,:));
        plot(x_1(4,:));
        plot(x2(4,:));
        legend('Q1=diag(0.1,5,0.1,5)','Q2=Q1*10','Q3=Q1*100') 
        title('State: Unsprung Mass Velocity')
        ylabel('dot(z_u) (m/s)')
        xlabel('simulation steps (0.01s/step)')
        hold off
%% R
        figure(2)
        subplot(2,2,1)
        hold on
        plot(x_000(1,:));
        plot(x4(1,:));
        plot(x5(1,:));
        plot(x6(1,:));
        title('State: Suspension Stroke')
        ylabel('z_s-z_u (m)')
        xlabel('simulation steps (0.01s/step)')
        hold off
        
        subplot(2,2,2)
        hold on
        plot(x_000(2,:));
        plot(x4(2,:));
        plot(x5(2,:));
        plot(x6(2,:));
        title('State: Sprung Mass Velocity')
        ylabel('dot(z_s) (m/s)')
        xlabel('simulation steps (0.01s/step)')
        hold off
        
        subplot(2,2,3)
        hold on
        plot(x_000(3,:));
        plot(x4(3,:));
        plot(x5(3,:));
        plot(x6(3,:));
        title('State: Tyre Deflection')
        ylabel('zt (m)')
        xlabel('simulation steps (0.01s/step)')
        hold off
        
        subplot(2,2,4)
        hold on
        plot(x_000(4,:));
        plot(x4(4,:));
        plot(x5(4,:));
        plot(x6(4,:));
        legend('R1=diag(0.0000001,0.0000001)','R2=R1*0.1','R3=R1*0.01','R4=R1*0.001') 
        title('State: Unsprung Mass Velocity')
        ylabel('dot(z_u) (m/s)')
        xlabel('simulation steps (0.01s/step)')
        hold off
        
        figure(4)
        hold on
        plot(u_rec_000(1,:))
        plot(u_rec4(1,:))
        plot(u_rec5(1,:))
        plot(u_rec6(1,:))
        hold off
        title('Input: Actuator Force (N)')
        legend('R1=diag(0.0000001,0.0000001)','R2=R1*0.1','R3=R1*0.01','R4=R1*0.001')
        ylabel('Actuator Force (N)')
        xlabel('simulation steps (0.01s/step)')
        
        figure(5)
        plot(u_rec_000(2,:))
        title("Bump Road Preview Disturbance in Regulation MPC")
        ylabel('z_s (m/s)')
        xlabel('simulation steps (0.01s/step)')
%% N
        figure(3)
        subplot(2,2,1)
        hold on
        plot(x12(1,:));
        plot(x11(1,:));
        plot(x10(1,:));
        plot(x8(1,:));
        plot(x9(1,:));
        title('State: Suspension Stroke')
        ylabel('z_s-z_u (m)')
        xlabel('simulation steps (0.01s/step)')
        hold off
        
        subplot(2,2,2)
        hold on
        plot(x12(2,:));
        plot(x11(2,:));
        plot(x10(2,:));
        plot(x8(2,:));
        plot(x9(2,:));
        title('State: Sprung Mass Velocity')
        ylabel('dot(z_s) (m/s)')
        xlabel('simulation steps (0.01s/step)')
        hold off
        
        subplot(2,2,3)
        hold on
        plot(x12(3,:));
        plot(x11(3,:));
        plot(x10(3,:));
        plot(x8(3,:));
        plot(x9(3,:));
        title('State: Tyre Deflection')
        ylabel('zt (m)')
        xlabel('simulation steps (0.01s/step)')
        hold off
        
        subplot(2,2,4)
        hold on
        plot(x12(4,:));
        plot(x11(4,:));
        plot(x10(4,:));
        plot(x8(4,:));
        plot(x9(4,:));
        legend('N=2','N=4','N=6','N=8','N=15') 
        title('State: Unsprung Mass Velocity')
        ylabel('dot(z_u) (m/s)')
        xlabel('simulation steps (0.01s/step)')
        hold off