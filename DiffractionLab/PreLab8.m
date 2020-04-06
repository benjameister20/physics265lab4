R=0.45e-3; %radius of the aperture
lambda=632.8e-9; % wavelength of the laser
Xs=[-0.5:0.001:0.5]*1e-3; % vector of distances to the center of the screen.
                          % the distance between points is 1 micron
d=20e-2; % distance from the aperture to the screen

% Airy disk pattern. I introduce the variable x here for two purposes.
% First, as you see, there are two places in line 15 where x is used. 
% If I replace x in line 15 back to its actual equation, right side of 
% line 14, then Matlab will calculate it twice. So, calculating x
% separately makes code faster.
% Second, it makes line 15 to be less crowded and easier to read.

x=2*pi*R/lambda*Xs./sqrt(Xs.^2+d^2);
I=(2*besselj(1,x)./x).^2;
% I made I0 to be equal to 1.

% Plotting the result
figure(1);
plot(Xs,I);