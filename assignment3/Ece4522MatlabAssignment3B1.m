% ./Ece4522/MatlabAssignment3/B1.m
% Edge Detection and Location in 1-D Signals
%
% By: Leomar Duran <https://github.com/lduran2/>
% When: 2021-11-03t20:53
% For: ECE 4522
% Version: 1.3.0
%
% CHANGELOG:
%     v1.3.0 - 2021-11-03t20:53
%         compressed the threshold signal into location signal and
%         plotted it
%
%     v1.2.0 - 2021-11-03t19:46
%         add the thresholding
%
%     v1.1.1 - 2021-11-03t17:34
%         formatted the plots
%
%     v1.1 - 2021-11-03t17:11
%         plotted x[n] and y[n]
%
%     v1.0 - 2021-11-03t16:49
%         finding the filter for backward difference system
%

clear;
%% Overview of a Backward Difference System

N = 159;    % number of samples
nn = 1:N;   % the samples
amp = 255;  % the amplitude
xx = amp*(rem(nn,30)>19);    % the signal x[n]

% A discrete-time backward difference system has the input-output
% relation given by
%   y[n] := x[n] - x[n - 1]
%   <=>         Y(z) := X(z) - X(z)z^-1 = X(z)(1 - z^-1)
%    => H(z) := Y(z)/X(z)               =     (1 - z^-1)/1
B = [1, -1]; % (1 - z^-1)
A = [1]; % 1

yy = filter(B, A, xx);  % digital filter with impulse response
                        % x[n] -> y[n]

%% thresholding the backward difference
T = amp * 0.6           % display threshold at 60% of the amplitude
dd = (abs(yy) >= T);    % the threshold signal

%% find the nonzero output locations
LL = find(dd)   % display the locations

%% plot the input and output signal
figure(1);

funcs = [ xx; yy; dd ];         % functions to graph
nfuncs = (size(funcs,1) + 1)    % display number of functions
titles = {
    '(a) input signal x[n]',
    '(b) output signal of backward difference system y[n]',
    '(c) thresholding'
    };

% limits
% display these for debugging
ymin = -amp
ymax = round((amp + 50), -2)
ymins = [ repmat(ymin,1,2), 0 ]
ymaxs = [ repmat(ymax,1,2), 1 ]

ysignalLabel = 'magnitude [<1>]';
ylabels = { ysignalLabel, ysignalLabel, 'is edge? [F, T]' };

for k = 1:(nfuncs - 1)
    subplot(nfuncs, 1, k);
    stem(nn, funcs(k,:));

    xlim([0, N]);
    ylim([ymins(k), ymaxs(k)]);

    title(titles{k});
    xlabel('time index [sample]');
    ylabel(ylabels{k});
end % for k = 1:nfuncs

%% the location signal is a little different to plot
subplot(nfuncs, 1, nfuncs);
lenLL = length(LL)          % display length of LL signal
stem(LL, ones(1,lenLL));    % plot LL to 1s

xlim([0, N]);   % match transpose axis of the signal above

title('(d) edge locations');
xlabel('edge location [sample]');
ylabel('1 if edge');

disp('Done.');