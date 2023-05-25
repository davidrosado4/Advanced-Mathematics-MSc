% David Rosado Rodríguez

%---------------------------------EXERCISE 1------------------------------------

% Read the wav file
[x,FS] = audioread("easy.wav");


% Play the audio
player = audioplayer (0.8*x, FS);
play (player);


% Add noise
y = x + normrnd(0, 0.2, length(x),1);


% Wait for the first audio file to finish playing
while isplaying(player)
    pause(0.8);
end

% Play the noisy audio
player = audioplayer (y, FS);
play (player);

% Remove noise of the audio

% Preservation ratio
r=0.05;

% Analysis filters
wa='db6';
ws = 'db6';

% No. of levels
J=5;

% Plot the fwt of the original audio
figure(1);
plot(fwt(x,wa,J));

% Remove the small coefficients in the wavalet space

% FWT
c_fwt = fwt(y,wa,J);

% Plot it
figure(2);
plot(c_fwt);

% Keep the largest
[cc_fwt,n]=largestr(c_fwt,r);

% Plot it
figure(3);
plot(cc_fwt);

% Go back the audio
r_fwt =ifwt(cc_fwt,ws,J,length(y));

% Wait for the noisy audio file to finish playing
while isplaying(player)
    pause(0.8);
end

% Play the new audio
player = audioplayer (r_fwt, FS);
play (player);



