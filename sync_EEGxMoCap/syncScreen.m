%% Screen modification
% Get the screen size
screenSize = get(0, 'ScreenSize');
width = screenSize(3); % Width of the screen
height = screenSize(4); % Height of the screen

% Create a black image matrix
blackImage = zeros(height, width);
whiteImage = ones(height,width);

%% Pause times
% Define the number of delay times 
numDelays = 10; 

% Define the range for the delay times (in seconds)
minDelay = 0.1; % Minimum delay time (in seconds)
maxDelay = 1; % Maximum delay time (in seconds)

% Generate random delay times array
delayTimes = (maxDelay - minDelay) .* rand(1, numDelays) + minDelay;

%% Displaying screens with delays
fig = figure(1);
b = 0; % if b = 0, screen should be black, b = 1 - screen is white
i = 1;

while i < length(delayTimes)
    pause(delayTimes(i))
    fig.WindowState = 'fullscreen';
    fig.Units = 'normalized';
    fig.OuterPosition = [0 0 1 1];
    if (b)
        imshow(whiteImage,'Border','tight');
    else 
        imshow(blackImage,'Border','tight');
    end

    screenCol(i) = b; 
    i = i + 1;
    b = ~b;
end

% setup EEG and camera software so it's ready to go
% start recording
% start screen script

% the screen code:
% show image at the start to indicate start of recording is coming - maybe
% a countdown?
% generate random number of delays and store in time variable
% alternate between black and white screens


