%function info=speckle2D(dataFile)
%speckle pattern 2d

%lambda=info.lambda;%wavelength [m]
clear;



% --------------- Load photon exit location and phase data ----------------
dataFile = 'exit-aperture-101.txt';
data = dlmread(dataFile);

lambda = 780e-7;
 
% Distance between medium and detector.
D = 5; % [cm]



%define camera
CCDGrid=zeros(100,100);
CCDdx=2.5e-3/size(CCDGrid,1);
CCDdy=2.5e-3/size(CCDGrid,2);

% The aperture of the medium (window in which photons leave)
% is defined as a 3x3cm square centered at x=5,y=5.  Therefore,
% it extends to 3->7cm.  In order to calculate the distance from
% each photon exit location to pixel on the CCD camera, we must
% know the location of the pixel in 2-D space.  It is assumed
% that the CCD camera, regardless of it's size is centered at the
% same location as the aperture.  That is, it's midpoint is at
% x=5,y=5.
start_x = 5 - (size(CCDGrid,1)/2*CCDdx);
start_y = 5 - (size(CCDGrid,2)/2*CCDdy);

% only grab a chunk of photons for testing.
data = data(1:100,:);

% for each pixel in the x-axis
for i = 1:size(CCDGrid, 1)
    x_pixel = start_x + (CCDdx * i);
    
    % for each pixel in the y-axis
    for j = 1:size(CCDGrid, 2)
        y_pixel = start_y + (CCDdy * j);
        
        % for each photon
        for k = 1:size(data, 1)
            
            % retrieve the weight of the photon.  
            weight = data(k,1);
            
            % transmission angle of the photon when it exiting
            % the medium.
            transmission_angle = cos(data(k, 2));
            
            % retrieve the stored path length in the medium.
            path_length = data(k,3);
            
            % retrieve the exit location of the photon.
            x = data(k,4);
            y = data(k,5);
            
            
            
            
            
            % add the distance from medium to the detector pixel.
            dist_to_pixel = sqrt(D^2 + (x_pixel - x)^2 + (y_pixel - y)^2);
            L = path_length + dist_to_pixel;
            
            
            
            CCDGrid(i,j) = (CCDGrid(i,j) + 1/dist_to_pixel*weight*exp(-(1i)*L*2*pi/lambda));
        end
        
    end
    y_pixel = start_y + CCDdy;
end

figure;
imagesc((abs(CCDGrid)))
colormap hot
drawnow
            


% function Norms=getNorms(stuff)
% tic
%     CCDdx=stuff.CCDdx;
%     CCD2Surf=stuff.CCD2Surf;
%     CCDGrid=stuff.CCDGrid;
%     surfaceGrid=stuff.surfaceGrid;
%     surfacedX=stuff.surfacedX;
%     lambda=stuff.lambda;
% 
% offsetx=((size(CCDGrid,1)-0.5)*CCDdx-(size(surfaceGrid,1)-0.5)*surfacedX)/1;
% offsety=((size(CCDGrid,2)-0.5)*CCDdx-(size(surfaceGrid,2)-0.5)*surfacedX)/1;
% norm1=zeros(size(surfaceGrid));
% %Norms{1:size(CCDGrid,1),1:size(CCDGrid,2)}=zeros(size(surfaceGrid));
%     for i=1:size(CCDGrid,1)
%         for j=1:size(CCDGrid,2)
%             y=[i*CCDdx j*CCDdx CCD2Surf];
%             
%             for k=1:size(surfaceGrid,1);
%                 for l=1:size(surfaceGrid,2);
%                     x=[k*surfacedX-offsetx l*surfacedX-offsety 0];
%                     norm1(k,l)=norm(y-x);
%                 end
%             end
%             Norms{i,j}=exp(1i*2*pi/lambda*norm1)./norm1;
%         end
%     end
%     toc
%     disp ('geometry matices created')
%  
% end