% Second Laboratory of Applied Harmonic Analysis
% David Rosado Rodríguez

% ------------------------------SOLUTION----------------------------------------

% Three audios will sound. The first one is the normal audio. The second one is 
% the noisy one and the last one is the processed audio

% Read the wav file
[x,FS] = audioread("easy.wav");

% Plot frequencies
figure;
freqz(x);

% Play the audio
player = audioplayer (0.8*x, FS);
play (player);

% Add noise 
t = (1:length(x))/FS;
y = 0.8*x + 0.1*sin(35000* t');

% Plot frequencies
figure;
freqz(y);

% Wait for the first audio file to finish playing
while isplaying(player)
    pause(0.8);
end

% Play the noisy audio
player = audioplayer (y, FS);
play (player);

% Remove noise using the filter h 
theta = 35000/FS;
b = zeros(1,3);
b(1) = 1; b(2) =-2*cos(theta); b(3) = 1;
a = zeros(1,3);
a(1) = 1;
yy = filter(b,a,y);

% Plot frequencies
figure;
freqz(yy);
% Wait for the second audio file to finish playing
while isplaying(player)
    pause(0.8);
end

% Play the noisy audio
player = audioplayer (yy ,FS);
play (player);

% Finally, let us compute and plot the magnitude of the filter’s frequency response
b(1) = 1; b(2) =-2*cos(pi/6); b(3) = 1;
figure;
freqz(b);





