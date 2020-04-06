%% Section 1. Variables and parameters
% laser wavelength
lambda = 632.8*10^-9; % m

% distance between aperture and screen
d = 20*10^-2; % m

% aperture radius
Aperture_Radius = (0.9/2)*10^-3; % m

% pixel width
pixel_width = 4.8*10-6; % m

% screen size
Screen_Size = 1000*pixel_width; % m

% aperture area
Aperture_Area = pi*(Aperture_Radius)^2; % m^2

% wavenumber
k = 2*pi/lambda;

fprintf("lambda = %d m", lambda);
fprintf("d = %d m", d);
fprintf("Aperture_Radius = %d m", Aperture_Radius);
fprintf("pixel_width = %d m", pixel_width);
fprintf("Screen_Size = %d m", Screen_Size);
fprintf("Aperture_Area = %d m^2", lambda);
fprintf("k = %d", k);

%% Section 2. Dividing aperture into area elements
% Calculate the initial value dxa for the spacing between aperture elements
d_a = Aperture_Radius;
d_s = 1/2*Screen_Size;
dxa = lambda/(2*(d_s + d_a)) * sqrt((d_s + d_a)^2 + d^2); %

% For large distances between the apperture and the target, the distance
% dxa might be comparable, or even larger than the aperture size. Here I
% check if the initial distance is large and limit it to 0.1 mm
if dxa > 0.1e-3
    dxa = 0.1e-3;
end

%If you want to make dxa smaller, then you can uncomment next line and use it.
%dxa=dxa/10;

N = ceil(Aperture_Radius/dxa); % Number of  elements that fits along the aperture radius;
Xa = -N*dxa:dxa:N*dxa; % Vector that contains the x coordinates of elements

[XaM, YaM] = meshgrid(Xa,Xa); % Matrix XaM contains X coordinates of all elements;
                              % Matrix YaM contains Y coordinates of all elements;

%% Section 3. Initializing the amplitude of the wave A(r_0) at the aperture
% Create a 2D array of amplitudes of the same size as arrays XaM and YaM of x and y coordinates
% All values are set to 0
A0 = zeros(size(XaM));

% Calculate the amplitude of the wave here:
P_a = 1; % W
Screen_Area = (1024*4.8*10^-6) * (1280*4.8*10^-6); % m^2
A_0 = sqrt(P_a/Screen_Area); % W/m^2

% 2D array that contains distances from corresponding elements to the
% center of the aperture;
Ra = sqrt(XaM.^2+YaM.^2);

% 2D logical array. Element of this array is equal to 1 (true) ff distance of the corresponding 
% aperture element to the center of the aperture is less than the aperture radius (element belongs to the aperture).
ind = (Ra<=Aperture_Radius);

% Setting amplitudes of all elements that are inside the aperture to the
% initial value
A0(ind) = A_0;

% Visualization of the amplitude array, just in case
figure(1);
imshow(A0./max(A0(:)));

%% Section 4. Dividing screen into area elements
% Calculate the initial value dxs for the spacing between screen elements elements
dxs = lambda/(2*(d_s+d_a))*sqrt((d_s+d_a)^2+d^2);

% Use what you have learned in Section 2 and create here a 1D array that contains x coordinated 
% of elements
Xs = -d_s:4.8*10-6:d_s;

%% Section 5. Constructing diffraction pattern
% This section is the main body of the simulation. In this section you have to:
% 1) Calculate the complex amplitude A(r) of the wave and its intensity I_wave(r) for every element of the array of screen elements Xs;
    % wavenumber
    k = 2*pi/lambda;

    K = atan(d./dxs);
    

% 2) calculate the Airy disk pattern I_Airy(r) for the array Xs, and compare it to I_wave(r);
    % angle of observation
    theta = 0:0.1:2*pi;

    % aperture radius
    a = 0.9/2*10^-3;

    % Bessel function of first kind of order one
    J_1 = besselj(1, k*a*sin(theta));

    % Initial intensity
    I_0 = 1;

    % Airy Disk intensity function
    I = I_0.*(2.*J_1./(k.*a.*sin(theta))).^2;

    figure(2);
    plot(theta, I);
    xlabel("Theta");
    ylabel("Intensity");

% 3) Create 2D arrays of x and y coordinates of the screen elements, analogues to XaM and YaM, see section 2 of the script;
    [XsM, YsM] = meshgrid(Xs,Xs); % Matrix XaM contains X coordinates of all elements;
                                  % Matrix YaM contains Y coordinates of all elements;

% 4) Calculate the intensity I_wave(r) for all elements on the screen, and plot it as a 2D image (analogue of the picture taken by camera);
    r_vals = sqrt(XsM.^2 + YsM.^2 + d.^2);
    syms r_a
    r = XsM;
    A = symsum(A_0.*K.*exp(1i.*k.*abs(r - r_a))./abs(r - r_a), r_a, XsM);
    I_wave = A.^2;
    figure(3);
    plot(r, I_wave);
    
% 5) Calculate the total power of the simulated diffraction pattern. It should be equal to the power of the wave at the aperture.
    S = [XsM, YsM];
    P_s = I_wave .* S;

%% Section 6. Extracting diffraction pattern from an image
IMAGE=imread('Image0.bmp');
figure(3);
imshow(IMAGE);

% Use ZOOM IN and DATA CURSOR tools of the figure window to find the
% center of the diffraction pattern.

% To be free of these restrictions we will convert uint8 type into double
% type:

IMAGE_dbl=double(IMAGE);

% Now you can calculate baground, rescale image data, extract cross section
% through the center of the diffraction pattern and compare it with the
% simulation



