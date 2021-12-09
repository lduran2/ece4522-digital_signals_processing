% ./Ece4522/MatlabAssignment3/B2OFF.m
% Bar Code Detection, Rotation and Decoding
%
% By: Leomar Duran <https://github.com/lduran2/>
% When: 2021-11-11t16:51
% For: ECE 4522
% Version: 3.0.0
%
% CHANGELOG:
%     v3.0.0 - 2021-11-11t16:51
%         attempt at correcting and decoding a rotated barcode image
%
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
IM_FILENAME = 'OFFv3.png';
rotimdata = imread(IM_FILENAME);    % data read from the image
[rotNrows, rotN] = size(rotimdata); % the dimensions of the image

% find the slope rows, assuming 100% vertical padding
islope_rows = round(rotNrows*[1 2]/3)  % display rows to scan rotation
slopes_n = rotimdata(islope_rows,:);   % slope rows

amp = double(max(max(abs(slopes_n))))   % display the amplitude
T = amp * 0.6                   % display threshold at 60% of the amplitude

%% Overview of a Backward Difference System
% A discrete-time backward difference system has the input-output
% relation given by
%   y[n] := x[n] - x[n - 1]
%   <=>         Y(z) := X(z) - X(z)z^-1 = X(z)(1 - z^-1)
%    => H(z) := Y(z)/X(z)               =     (1 - z^-1)/1
B = [1, -1]; % (1 - z^-1)
A = [1]; % 1

% apply the system to the two slope rows first
filter_slopes_n = @(k) filter(B, A, slopes_n(k,:));
% pass each row individually, resulting in horizontal table
yslopes_n_cells = arrayfun(filter_slopes_n, [1 2], 'UniformOutput', false);
% flip table and convert to a matrix
yslopes_n = cell2mat(yslopes_n_cells');

%% threshold the backward difference of the slope rows
dslopes_n = (abs(yslopes_n) >= T);  % the threshold of the slope rows

%% find the first edge of each slope row
% ignore the first set of points
skip_to = 4
find_slopes_skip = @(k) find(dslopes_n(k,:),skip_to);
Lslopes_skip_cells = arrayfun(find_slopes_skip,[1:2],'UniformOutput',false)
Lslopes_skip_1_pad = cell2mat(Lslopes_skip_cells')
Lslopes_skip = Lslopes_skip_1_pad(:,skip_to)

% calculate the angle of rotation
slope_d = [Lslopes_skip'; islope_rows] * [-1; 1]    % dx,dy
slope_atan = atan2(slope_d(2), slope_d(1))          % r [rad]
slope_deg = slope_atan*360/(2*pi)                   % r [deg]

unrot_atan = (slope_atan - pi/2)    % corresponding unrotate angle
                                    % negative complement

% rotation function
RT = { @(T) cos(T), @(T) sin(T); @(T) -sin(T), @(T) cos(T) };
Runrot = cellfun(@(f) f(unrot_atan), RT)    % rotate with the unrotate
                                            % angle

% find the original corners
pcorners = [ 0, rotN ];        % original corners x-coord
qcorners = [ 0, rotNrows ];    % original corners y-coord
% find the cartesian product of the corner coords
[xcorners, ycorners] = meshgrid(pcorners, qcorners);
corners = [xcorners(:), ycorners(:)]
% indices from 1 to the # of corners
icorners = 1:size(corners, 1);

Rxcorners = @(k) Runrot*corners(k,:)';      % corner rotating function
% apply to each corner, the result will be a table of 2x1 matrices
unrotcorners = arrayfun(Rxcorners, icorners, 'UniformOutput', false);
% find the maximum of each coordinate, this is the size of the image
unrotsize = flip(round(max(cell2mat(unrotcorners)')+1));
nrows = unrotsize(1)
N = unrotsize(2)

% white out an array
imdata = ones(nrows, N) * 255;

% copy over the matrix
% note that matrix is irow (y) x icol (x)
for irow = 1:rotNrows
    for icol = 1:rotN
        % find the unrotated coordinates
        coords = round(abs(Runrot*([icol; irow] - 1))) + 1';
        % if inside the range of index
        if (all(coords >= 1))
            % map the pixel
            imdata(coords(2), coords(1)) = rotimdata(irow, icol);
        end % if (all(coords >= 1))
    end % for icol = 1:N
end % for irow = 1:nrows

% display the unrotated image
figure(1)
imagesc(imdata);
colormap([ zeros(1,3); ones(1,3) ]);

%% now find the input row
nn = 1:N;                       % the samples
ixxrow = round(nrows/2)         % display index of row to use for xx
xn = imdata(ixxrow,:);          % the input data (middle row of image
                                % data)

%% backward difference on the input row
yn = filter(B, A, xn);  % digital filter with impulse response
                        % x[n] -> y[n]

%% thresholding the backward difference of the input row
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
for irow = ee
    if (irow <= midBarMin)
        loLe = 1
        upLe = nbars
    elseif (irow >= midBarMax)
        loLe = loLeMax
        upLe = lenLe
    else
        loLe = (irow - midBarMin + 1)
        hiLL = (loLe + nbars - 1)
    end
    % to get the basic width,
    % divide the first (nbars) edges by (nwidths)
    w1 = ((Le(upLe) - Le(loLe)) / ntotalWidths) % display the basic width
    NWIDTH(irow) = round(NWIDTH(irow) / w1);          % get number of widths
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
figure(2);

% limits
% display these for debugging
ymin = -amp
ymax = round((amp + 50), -2)
ymins = [ repmat(ymin,1,2), 0 ]
ymaxs = [ repmat(ymax,1,2), 1 ]

for irow = 1:(nnnfuncs - 1)
    subplot(nnnfuncs, 1, irow);
    stem(nn, nnfuncs(irow,:));

    xlim([0, N]);
    ylim([ymins(irow), ymaxs(irow)]);

    title(nntitles{irow});
    xlabel('time index [sample]');
    ylabel(nnylabels{irow});
end % for k = 1:(nfuncs - 1)

%% the location signal is a little different to plot
subplot(nnnfuncs, 1, nnnfuncs);
stem(Le, ones(1,lenLe));    % plot Le to 1s

xlim([0, N]);   % match transpose axis of the signal above

title(['(d)' eetitles{1}]);
xlabel(eeylabels{1});
ylabel('1 if edge');

%% plot the edge index domain functions
figure(3)

for irow = 1:(neefuncs - 1)
    subplot(neefuncs, 1, irow);
    stem(ee, eefuncs(irow,:));
    title(eetitles{irow});
    xlabel(eexlabel);
    ylabel(eeylabels{irow});
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
