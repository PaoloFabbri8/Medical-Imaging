%% Medical Imaging 2026 Exercise 2: MRI & NIRS
%
% This demo illustrates the features of the k-space.
% No coding needed, just follow the instructions in the exercise sheet
%
% Paavo Hietala 2026

clear; clc; close all;
    
%% 2. Generate a phantom and the k-space representation

N = 256; % image size
center = N/2;

% Add Gaussian noise to simulate scanner noise

noise_std = 0.02;

try
    P = phantom(N);  % Phantom; our "body" to be imaged
catch
    warning('No Image Processing Toolbox -- Loading the .MAT version');
    load('phantom_256.mat');
end

phantom_img = P + noise_std*randn(N,N);

% Compute full k-space
kspace = fftshift(fft2(phantom_img));

figure;
subplot(1,2,1);
imshow(phantom_img, []);
title('Original Phantom Image');
xlabel('x'); ylabel('y');

subplot(1,2,2);
imagesc(log(abs(kspace)));
colormap gray; axis image;
title('Log-Magnitude of k-space');
xlabel('kx'); ylabel('ky');
    
%% 3. Plot a couple of k-space points
%
% Note! The k-space values are complex numbers, here we are looking at
% the corresponding cosine wave.

N = 256;
[x, y] = meshgrid(0:N-1, 0:N-1);

% Shift coordinates to center
x_c = x - N/2;
y_c = y - N/2;

% k-space points
points = [10, 10; -10, -10; 38, 7]; 
colors = {'r','g','b'};

figure;

% Top-left: k-space plot with marked points
subplot(2,2,1); hold on;
imagesc(log(abs(kspace)+1)); colormap gray; axis image;
title('Full k-space (log-magnitude)');
xlabel('kx'); ylabel('ky');

for idx = 1:size(points,1)
    kx = points(idx,1) + N/2 + 1;
    ky = points(idx,2) + N/2 + 1;
    plot(kx, ky, 'o', ...
        'MarkerEdgeColor', colors{idx}, ...
        'MarkerFaceColor', colors{idx}, ...
        'MarkerSize', 8);
end
legend({'Point 1','Point 2','Point 3'});

% Remaining subplots: individual waves
for idx = 1:size(points,1)
    kx = points(idx,1);
    ky = points(idx,2);

    wave = cos(2*pi*(kx*x_c/N + ky*y_c/N));
    
    subplot(2,2,idx+1);
    imagesc(wave); colormap gray; axis image;
    title(sprintf('Wave from k-space point (%d,%d)', kx, ky));

    set(gca,'YDir','normal');
end
%% 4. Remove the center of k-space
%
% Feel free to play around with radius_center and see what happens!

k_mask_center = kspace;
radius_center = 10; % pixels to block

[x,y] = meshgrid(1:N,1:N);

mask_center = sqrt((x-center).^2 + (y-center).^2) <= radius_center;
k_mask_center(mask_center) = 0;

img_mask_center = abs(ifft2(ifftshift(k_mask_center)));

figure;
subplot(1,2,1);
imagesc(log(abs(k_mask_center)));
colormap gray; axis image;
xlabel('kx'); ylabel('ky');
subplot(1,2,2);
imshow(img_mask_center, []);
title('Image after blocking center of k-space');
xlabel('x'); ylabel('y');
    
%% 5. Remove the edges of k-space
%
% Feel free to play around with radius_edges and see what happens!

k_mask_edges = kspace;
radius_edges = 90; % pixels to keep

mask_edges = sqrt((x-center).^2 + (y-center).^2) > radius_edges;
k_mask_edges(mask_edges) = 0;

img_mask_edges = abs(ifft2(ifftshift(k_mask_edges)));

figure;
subplot(1,2,1);
imagesc(log(abs(k_mask_edges)));
colormap gray; axis image;
xlabel('kx'); ylabel('ky');
subplot(1,2,2);
imshow(img_mask_edges, []);
title('Image after blocking edges of k-space');
xlabel('x'); ylabel('y');
    
%% 6. Half-Fourier (upper half of the k-space)

% Zero the bottom half of the k-space
k_half = kspace;
k_half(center+1:N,:) = 0;

% Plot the image from half-filled k-space
img_half = abs(ifft2(ifftshift(k_half)));

figure;
subplot(2,2,1);
imagesc(log(abs(k_half)));
colormap gray; axis image;
xlabel('kx'); ylabel('ky');
subplot(2,2,2);
imshow(img_half, []);
title('Image from half k-space (top half removed)');
xlabel('x'); ylabel('y');

% Compute the bottom half of the k-space based on acquired data
k_unshifted = ifftshift(kspace);

k_half = zeros(size(kspace));
k_half(1:center,:) = k_unshifted(1:center,:);

rows_lower = center+1:N;
rows_upper = N - rows_lower + 2;  % mirrored row indices

cols = 1:N;
cols_mirror = N - cols + 2;       % mirrored column indices
cols_mirror(1) = 1;               % special case j=1

% Fill lower half using conjugate symmetry
k_half(rows_lower, :) = conj(k_half(rows_upper, cols_mirror));

k_half = ifftshift(k_half);

% Reconstruct image from mirrored k-space
img_half = abs(ifft2(k_half));

subplot(2,2,3);
imagesc(log(abs(k_half)));
colormap gray; axis image;
xlabel('kx'); ylabel('ky');
subplot(2,2,4);
imshow(img_half, []);
title('Corrected Half-Fourier');
xlabel('x'); ylabel('y');

%% Bonus: Animated reconstruction from k-space points

N = 128; % smaller for faster animation
phantom_img = phantom(N);

% Compute k-space
kspace = fftshift(fft2(phantom_img));

% Prepare grid
[x, y] = meshgrid(0:N-1, 0:N-1);

% Initialize sum image
img_sum = zeros(N);

% Get coordinates of k-space points in linear order (flattened)
[kx_list, ky_list] = meshgrid(1:N, 1:N);

% Flatten k-space indices and corresponding values
kx_list = kx_list(:);
ky_list = ky_list(:);
kspace_vals = kspace(:);

% Sort by magnitude to add strongest contributors first
[~, idx_sort] = sort(abs(kspace_vals), 'descend');
kx_list = kx_list(idx_sort);
ky_list = ky_list(idx_sort);
kspace_vals = kspace_vals(idx_sort);

figure;
h = imagesc(img_sum); colormap gray; axis image; axis off;
title('Fourier Synthesis: Adding k-space points');
drawnow;

% Animation loop (add top 500 points for speed)
for n = 1:500
    kx = kx_list(n) - N/2 - 1; % center at 0
    ky = ky_list(n) - N/2 - 1;
    A = kspace_vals(n);
    
    % Generate 2D sinusoid
    wave = real(A * exp(1i*2*pi*(kx*x/N + ky*y/N)));
    
    % Add to cumulative image
    img_sum = img_sum + wave;
    
    % Update plot
    set(h,'CData',img_sum);
    title(sprintf('Added %d k-space points', n));
    drawnow;
    
    % Pause briefly to visualize
    pause(0.01);
end