function [dot_zr] = bumproad(t,ts,t0,V,L,A)

% Based on Agostinacchio, M., D. Ciampa, and S. Olita. "The vibrations 
% induced by surface irregularities in road pavements�a Matlab� approach." 
% European Transport Research Review 6.3 (2014): 267-275.

%% Variables
% t:   The time length of the bump road;
% ts:  The sampling time of the bump road;
% t0:  The starting time of the bump;
% V:   The velocity of the car;
% L:   The length of the bump;
% A:   The height of the bump; 
dot_zr=zeros(1,t/ts);
for i=ts:ts:t
    if i>= t0 && i<(L/V)
        dot_zr(round(i/ts))=A/2*(1-cos(2*pi*V*(i-t0)/L));
    end
end