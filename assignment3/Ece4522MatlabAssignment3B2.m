% ./Ece4522/MatlabAssignment3/B2.m
% Bar Code Detection and Decoding
%
% By: Leomar Duran <https://github.com/lduran2/>
% When: 2021-11-04t03:12
% For: ECE 4522
% Version: 2.3.0
%
% CHANGELOG:
%     v2.3.0 - 2021-11-04t03:12
%         calculated the encoded message, properly decoded it
%
%     v2.2.0 - 2021-11-04t00:16
%         calculated basic width
%
%     v2.1.0 - 2021-11-03t23:35
%         found difference of locations, plotted it
%         cleaned up a bit
%
%     v2.0.0 - 2021-11-03t21:35
%         adapted to read from the image file
%
%     v1.3.0 - 2021-11-03t20:53
%         compressed the threshold signal into location signal,
%         plotted it
%
%     v1.2.0 - 2021-11-03t19:46
%         add the thresholding
%
%     v1.1.1 - 2021-11-03t17:34
%         formatted the plots
%
%     v1.1 - 2021-11-03t17:11
%         plotted x[n], y[n]
%
%     v1.0 - 2021-11-03t16:49
%         finding the filter for backward difference system
%

clear;
%% reading the image
IM_FILENAME = 'HP110v3.png';
imdata = imread(IM_FILENAME);   % data read from the image
imsize = size(imdata);          % the size of the image
ixxrow = round(imsize(1)/2)     % display index of row to use for xx
xn = imdata(ixxrow,:);          % the input data (middle row of image
                                % data)
amp = double(max(abs(xn)))      % display the amplitude
N = imsize(2)                   % display number of samples

%% Overview of a Backward Difference System
nn = 1:N;                   % the samples

% A discrete-time backward difference system has the input-output
% relation given by
%   y[n] := x[n] - x[n - 1]
%   <=>         Y(z) := X(z) - X(z)z^-1 = X(z)(1 - z^-1)
%    => H(z) := Y(z)/X(z)               =     (1 - z^-1)/1
B = [1, -1]; % (1 - z^-1)
A = [1]; % 1

yn = filter(B, A, xn);  % digital filter with impulse response
                        % x[n] -> y[n]

%% thresholding the backward difference
T = amp * 0.6           % display threshold at 60% of the amplitude
dn = (abs(yn) >= T);    % the threshold signal

%% find the nonzero output locations
Le = find(dn)       % display the location signal
lenLe = length(Le)  % display length of LL signal
ee = 1:lenLe;       % edge domain

%% backward difference on location signal
DELTAe = filter(B, A, Le)    % display the difference of location signal

%% find the basic width and # basic widths per bar

ndigits = 12;       % number of digits per signal
digit_bars = 4;     % number of bars per digit
digit_width = 7;    % total width of each digit

%   the number of (end + separator) bars
% =        2 ends * 3 widths/end
% | + 1 separator * 5 widths/separator
n_delimiter_bars = (2*3) + (1*5)    % display

nbars = ((ndigits * digit_bars) + n_delimiter_bars) % display total
                                                    % number of bars
midBarMin = round(((nbars - 1)/2) + 1)  % mid point to nbars
midBarMax = (nbars - midBarMin + 1)     % mid point to last nbars
loLeMax = (lenLe - nbars + 1)           % maximum lower bound of LL

%   total number of widths
% =     12 digits * 7 widths/digit
% | +      2 ends * 3 widths/end
% | + 1 separator * 5 widths/separator
ntotalWidths = ((ndigits * digit_width) + n_delimiter_bars)  % display

% calculate basic width for all values of edge index
NWIDTH = DELTAe; % copy delta signal to a number of signals symbol
for k = ee
    if (k <= midBarMin)
        loLe = 1
        upLe = nbars
    elseif (k >= midBarMax)
        loLe = loLeMax
        upLe = lenLe
    else
        loLe = (k - midBarMin + 1)
        hiLL = (loLe + nbars - 1)
    end
    % to get the basic width,
    % divide the first (nbars) edges by (nwidths)
    w1 = ((Le(upLe) - Le(loLe)) / ntotalWidths) % display the basic width
    NWIDTH(k) = round(NWIDTH(k) / w1);          % get number of widths
                                                % per DELTA
end % for k = ee

NWIDTH  % display the width signal

%% seek the beginning of the message
for iMsgStart = 1:(lenLe - 2)
    if (isequal(NWIDTH(iMsgStart:(iMsgStart + 2)), [ 1 1 1 ]))
        break
    end
end % for iMsgStart = ee

iMsg = (iMsgStart:(iMsgStart + nbars - 1))  % display message indices
MESSAGE = NWIDTH(iMsg)                      % display the message

%% functions and their labels

% time index domain

nnfuncs = [ xn; yn; dn ];           % functions to graph
nnnfuncs = (size(nnfuncs,1) + 1)    % display number of functions
nntitles = {
    '(a) input signal x[n]',
    '(b) output signal of backward difference system y[n]',
    '(c) thresholding, d[n]'
    };
ysignalLabel = 'magnitude [<1>]';
nnylabels = { ysignalLabel, ysignalLabel, 'is edge? [F, T]' };

% edge index domain

eefuncs = [ Le; DELTAe; NWIDTH ];   % functions to graph
neefuncs = (size(eefuncs,1) + 1)    % display number of functions
eetitles = { '(a) location signal l[n]', ...
    '(b) location difference signal \Delta[n]', '(c) number of w_1' };
eexlabel = 'edge index [edge]';
eeylabels = { 'edge location [sample]', 'bar width [sample]' , 'width [w_1]' };

%% plot the time index domain functions
figure(1);

% limits
% display these for debugging
ymin = -amp
ymax = round((amp + 50), -2)
ymins = [ repmat(ymin,1,2), 0 ]
ymaxs = [ repmat(ymax,1,2), 1 ]

for k = 1:(nnnfuncs - 1)
    subplot(nnnfuncs, 1, k);
    stem(nn, nnfuncs(k,:));

    xlim([0, N]);
    ylim([ymins(k), ymaxs(k)]);

    title(nntitles{k});
    xlabel('time index [sample]');
    ylabel(nnylabels{k});
end % for k = 1:(nfuncs - 1)

%% the location signal is a little different to plot
subplot(nnnfuncs, 1, nnnfuncs);
stem(Le, ones(1,lenLe));    % plot Le to 1s

xlim([0, N]);   % match transpose axis of the signal above

title(['(d)' eetitles{1}]);
xlabel(eeylabels{1});
ylabel('1 if edge');

%% plot the edge index domain functions
figure(2)

for k = 1:(neefuncs - 1)
    subplot(neefuncs, 1, k);
    stem(ee, eefuncs(k,:));
    title(eetitles{k});
    xlabel(eexlabel);
    ylabel(eeylabels{k});
end % for k = 1:(nfuncs - 1)

%% the message signal is also a little different to plot
subplot(neefuncs, 1, neefuncs);
stem(iMsg, MESSAGE);    % plot indices of message to message

xlim([0, round(lenLe, -1)]);   % match transpose axis of the signal above

title('(d) encoded message');
xlabel(eexlabel);
ylabel('encoded symbol');

%% decode the message
decodeMessage = decodeUPC(MESSAGE)  % decode

disp('Done.');
