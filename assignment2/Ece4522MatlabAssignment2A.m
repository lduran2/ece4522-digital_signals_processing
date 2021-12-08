% The Template from
% Part A. Overview of an FIR System

clear;

M = 99;
nn = 0:M; %<--discrete-time indices
xx = cos( 0.08*pi*nn ); %<--Input signal
bb = [1/3 1/3 1/3]; %<--Filter impulse response
yy = conv(bb, xx); %<--Compute the output

plot(nn, xx);
hold on
plot(nn, yy(nn+1));
