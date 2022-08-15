function [dot_zr] = roadGenerator2(roadClass,V,t,ts)

% Based on Agostinacchio, M., D. Ciampa, and S. Olita. "The vibrations 
% induced by surface irregularities in road pavements�a Matlab� approach." 
% European Transport Research Review 6.3 (2014): 267-275.

%% Variables
k = roadClass; % Values For ISO Road A-B Roughness Classification, from 3 to 9
number=t/ts; % Length of the Gaussian noise data
G0_library=[5;32;128;512;2048;8192;32768;131072]; % G_0 library from ISO 8608
G_0=G0_library(k,:)*10^(-6); % Corresponding G_0 value for road class
mu=0;
sigma=1;
noise=sigma*randn(1,number)+mu;  %whie Gaussian noise with mean=0 | covariance=1;

%% Road signal generaiton
dot_zr=2*pi*sqrt(G_0*V)*noise;


%% 
% % N  = 2500; %  Number of data points
% % L  = 250;  % Length Of Road Profile (m)
% B  = L/N ; % Sampling Interval (m)
% dn = 1/L;  % Frequency Band
% n0 = 0.1;        % Spatial Frequency (cycles/m)
% n  = dn : dn : N*dn; % Spatial Frequency Band
% 
% %% Road signal generation
% phi =  2*pi*rand(size(n)); % Random Phase Angle
% Amp1 = sqrt(dn)*(2^k)*(1e-3)*(n0./n); % Amplitude for Road  Class A-B
% x = 0:B:L-B; % Abscissa Variable from 0 to L
% hx = zeros(size(x));
% for i=1:length(x)
%     hx(i) = sum(Amp1.*cos(2*pi*n*x(i)+ phi));
% end
% figure
% plot(x,hx),xlabel('x (m)'),ylabel('h(x) (m)')
% 
% %% Road signal analysis
% [q , C] = psd_1D(hx, B, 'x');  % B is Sampling Interval (m); for the case that I have explained above it will be 250/45000 = 5.55e-3
% lambda = (2*pi) ./ q; % wavelengths
% f = q / (2*pi); % spatial frequency
% PSD = 2 * pi * C; % power spectrum
% figure
% loglog(f,PSD),xlabel('cycles/m'),ylabel('m^3')
% 
% xroad = x;
% hroad = hx;
end