% ./Ece4522/MatlabAssignment2/C.m
% FIR System
%
% By: Leomar Duran <https://github.com/lduran2/>
% When: 2021-10-14t04:29
% For: ECE 4522
% Version: 1.0
%
% CHANGELOG:
%     v1.0 - 2021-10-14t04:29
%         calculates the total system impulse response
%
%     v0.0 - 2021-10-13t22:54
%         template from part A

clear;
% parameters
q = 0.9;
r = 0.9;
M = 22;

nn = 0:M;       % time indices
d = [1];        % impulse signal
h1 = [1, -q];   % impulse response for w[.]
h2 = r.^nn;     % impulse response for y[.]

e = conv(h1, d);
f = conv(h2, e);

figure(1)
title("C. Cascading Two Systems")

stem(nn, f(nn+1));
xlabel('time index [sample]');
ylabel('magnitude <1>');

f(nn+1)
