% ./Ece4522/MatlabAssignment2/B1.m
% FIR System
%
% By: Leomar Duran <https://github.com/lduran2/>
% When: 2021-10-14t07:15
% For: ECE 4522
% Version: 1.4.1
%
% CHANGELOG:
%     v1.4.2 - 2021-12-08t19:08
%         playing audio subsequently, constants, commented to B.4
%             parameters
%
%     v1.4.1 - 2021-10-14t07:15
%         even vertical axes for figure 2
%
%     v1.4.0 - 2021-10-14t03:45
%         applied system, tested and write output file
%
%     v1.3.0 - 2021-10-14t03:17
%         implemented echo system impulse response, read input audio
%
%     v1.2.0 - 2021-10-14t01:31
%         plotted the error between x[n] and y[n]
%
%     v1.1.0 - 2021-10-14t01:15
%         created and plotted restoration system
%
%     v1.0.1 - 2021-10-13t23:36
%         moved plots to subplots, stem plot, labeled and ranged
%
%     v1.0 - 2021-10-13t23:12
%         defining the impulse response h[n]
%
%     v0.0 - 2021-10-13t22:54
%         template from part A
clear;

% constants
SHOW_PLOTS = true;  % whether to show the plots (for convenience)
PLAY_SOUNDS = 0b11; % how to play sounds
                    % (bit 0 = play input signal)
                    % (bit 1 = play output signal)

%% part B.1
M = 100;
nn = 0:M; %<--discrete-time indices
x = 256*(rem(nn,50)<10); %<--Input signal
h = [1, -0.9]; %<--Filter impulse response
w = conv(h, x); %<--Compute the output

if (SHOW_PLOTS)
    % Plot both the input and output signals on the same figure, using
    % `subplot` command.
    figure(1);
    title('B.1 FIR System');

    % the input signal plot
    subplot(2,1,1);
    stem(nn, x(nn+1));
    subtitle('Input signal, x[n]');
    xlim([0,75]);
    xlabel('time index [sample]');
    ylabel('magnitude <1>');

    % the output signal plot
    subplot(2,1,2);
    stem(nn, w(nn+1));
    subtitle('Output signal, w[n]');
    xlim([0,75]);
    xlabel('time index [sample]');
    ylabel('magnitude <1>');
end % if (showPlots)

%% part B.2 Restoration System
r = 0.9;
Mi = 22;
hi = r.^(0:Mi);

% Process the signal w[n] from Section B.1. with the system described
% by eq. (4) to obtain the output signal y[n].
y = conv(hi, w);

e = abs(y(nn+1) - x(nn+1)); % the error (difference) between y[n] and x[n]

y((23:32)+1) = y((23:32)+1);

if (SHOW_PLOTS)
    % Plot both the input and output signals on the same figure, using
    % `subplot` command.
    % Put the stem plots in side-by-side figures for comparison-use a
    % two-panel subplot.
    figure(2);
    title('B.2 Restoration System');

    % the input signal plot
    subplot(1,2,1);
    stem(nn, w(nn+1));
    subtitle('Input signal, w[n]');
    xlim([0,75]);
    ylim([-300,300]);
    xlabel('time index [sample]');
    ylabel('magnitude <1>');

    % the output signal plot
    subplot(1,2,2);
    stem(nn, y(nn+1));
    subtitle('Output signal, y[n]');
    xlim([0,75]);
    ylim([-300,300]);
    xlabel('time index [sample]');
    ylabel('magnitude <1>');

    % the error (difference) plot
    figure(3);
    title('B.2 Restoration System error (difference) plot, e[n]');
    stem(nn, e);
    xlim([0,50]);
    xlabel('time index [sample]');
    ylabel('magnitude <1>');
end % if (showPlots)

%% part B.3
% Evaluate the worst-case error by doing the following: use MATLAB's
% `max` function in the range n in [0, 50]
worstCaseErr = max(e((0:50)+1))

%% part B.4 An Echo System

% parameters
clear -regex [a-z]; % clear any variables, keeping constants
                    % (constant names have no lowercase letters)
%
r = 0.9;        % gain of the echo
P = 0.2*8000;   %[samples] echo delay = P [t] * sample requency FS
t_pad = 2;  %[s] extra wait between playing signals

% Implement the echo system in eq. (5) with the values of r and P
% determined in part (a).
h = zeros(1,P + 1);
h(1) = 1;
h(P + 1) = r;

% read and test the audio file
[Y,FS] = audioread('speech8k.wav');
if (bitand(PLAY_SOUNDS, 1))
    soundsc(Y);
end

% apply and test the echo system
Z = conv(h,Y')';
if (bitand(PLAY_SOUNDS, 2))
    % if also playing the input signal
    if (bitand(PLAY_SOUNDS, 1))
        % then wait for the audio to finish first
        duration_Y = (size(Y,1)/FS);
        pause(duration_Y + t_pad);
    end % if (bitand(PLAY_SOUNDS, 1))
    soundsc(Z);
end % if (bitand(PLAY_SOUNDS, 2))

% Save your result as a .wav file.
audiowrite('speech8k-echo.wav',Z,FS);

disp("Done.");
