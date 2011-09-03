%function info=speckle2D(dataFile)
% This script forms a speckle pattern on the face of a simulated CCD camera by
% interfering photons that exit from the Monte Carlo simulation.

%lambda=info.lambda;%wavelength [m]
clear;

% ------------------- Settings -------------------------------------
% Boolean values to decide if various graphs are formed.
createMovie = 1;
displayExitPhotons = 0;

% Wavelength of the photons.
lambda = 532e-7;

% Distance between medium and detector.
D = 40; % [cm]

% Start time of the simulation to make the speckle pattern from.
% That is, if the Monte Carlo + K-Wave simulation (i.e. the
% acousto-optic (AO) simulation) ran for a total time 'T', 'dt' is the time
% steps that occurred in the propagation of ultrasound.  Therefore we can
% take a portion of the whole AO simulation to make a speckle pattern, or
% use the entire photon exit data.
start_time = 200; % 'dt' starting at this time.
end_time = 200;  % end time that we want to look at.

% The acceptance angle of photons leaving the medium.
acceptance_angle = 0.90;



% The figure for the speckle pattern.
speckleFigure = figure;






%define the camera.
CCDGrid=zeros(200,200);
CCDdx=5.5e-2/size(CCDGrid,1);
CCDdy=5.5e-2/size(CCDGrid,2);

% The aperture of the medium (window in which photons leave)
% is defined in the Monte Carlo (MC) simulation.  Therefore,
% reference the simulation settings, specifically the detector size.
% In order to calculate the distance from
% each photon exit location to pixel on the CCD camera, we must
% know the location of the pixel in 2-D space.  It is assumed
% that the CCD camera, regardless of it's size is centered at the
% same location as the aperture.  That is, it's midpoint is at
% the same center coordinates as the exit aperture in the MC simulation.
%
center.x = 1.0; % Center of the CCD (which should be the center of exit aperture).
center.y = 1.0; % Note: in cm
start_x = center.x - (size(CCDGrid,1)/2*CCDdx);
start_y = center.y - (size(CCDGrid,2)/2*CCDdy);




% ------------------------- PLANE ---------------------------
% For drawing the plane made by the CCD and visualizing all
% photons that hit it.
%
% Define all the vertices on the plane (i.e. CCD corners)
if (displayExitPhotons)
    x = [start_x, (center.x+(center.x-start_x)), (center.x+(center.x-start_x)), start_x];
    y = [(center.y+(center.y-start_y)), (center.y+(center.y-start_y)), start_y, start_y];
    z = [D, D, D, D];
    color = [10, 10, 10, 10];
    rayFigure = figure;
    patch(x,y,z,color);
    view(3);
    % Define the exit aperture.
    t = 0:pi/360:2*pi;
    radius = 0.5;
    xy_plane = zeros(size(t)); xy_plane(:) = 2; % The xy-plane that the exit aperture resides on.
    color = zeros(size(t)); color(:) = 15; % Define the color for the exit aperture plane.
    hold on;
    patch(center.x + (cos(t)*radius), center.y + (sin(t)*radius), xy_plane, color);
    hold off;
end




% To create a movie that shows the fluctuating speckle patter (i.e.
% modulation) we loop over all exit data files, which contain information
% about the photons for a snapshot in time of the ultrasound.
if (createMovie)
    aviobj = avifile ( 'speckle-modulation2', 'compression', 'None', 'fps', 1);
end

for dt=start_time:end_time
    
    % Load the exit data from the AO simulation.
    dataFile = ['exit-aperture-', num2str(dt), '.txt'];
    data = dlmread(dataFile);
    
    display(sprintf('Time step (dt) = %i', dt));
    
    % -------------- Exit photon figure ----------------
    if (displayExitPhotons)
        % for each photon
        for photon = 1:size(data, 1)
            
            % retrieve the weight of the photon.
            weight = data(photon,1);
            
            % transmission angles of the photon when it exiting
            % the medium.
            dirx = data(photon, 2);
            diry = data(photon, 3);
            dirz = data(photon, 4);
            
            % Only plot photons that exit at a certain angle.
            if (dirz < acceptance_angle)
                continue;
            end
            
            % retrieve the stored path length in the medium.
            path_length = data(photon,5);
            
            % retrieve the exit location of the photon.
            x = data(photon,6);
            y = data(photon,7);
            z = data(photon,8);
            
            
            
            % Check for intersection of this photon with the CCD camera.
            %-----------------------------------------------------------
            % Exit location of photon.
            s0 = [x, y, z];
            % The normal to the CCD.
            %
            n = [0,0,1];
            
            % Point on the CCD.
            %
            p0 = [center.x, center.y, D];
            
            % Create the ray to cast from the exit of the photon out into
            % space.  If this ray intersects the CCD, then 'smear' it on
            % the face of the camera as if it was a wave.
            % NOTE:
            % - To extend the ray some portion beyond the plane made by the
            %   CCD, we give the 't' value (parametric form of line) z+D.
            %
            s1 = [x + dirx*(z+D), y + diry*(z+D), z + dirz*(z+D)];
            ray = [s0(1), s0(2), s0(3); ...
                s1(1), s1(2), s1(3)];
            
            figure(rayFigure);
            hold on;
            plot3(ray(:,1), ray(:,2), ray(:,3));
            hold off;
            
            % From the line segment (i.e. the ray 'r') we find the 't'
            % value that gives the intersection point on the CCD.
            %
            temp = p0 - s0;
            ti = dot(n, temp)/dot(n, s1);
            if (ti >= 0 & ti <= 1)
                xi = x + dirx*ti;
                yi = y + diry*ti;
                zi = D;
            end % if()
            
        end
    end
    
    % only grab a chunk of photons for testing.
    %data = data(1:500,:);
    
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
                
                % transmission angles of the photon when it exiting
                % the medium.
                dirx = data(k, 2);
                diry = data(k, 3);
                dirz = data(k, 4);
                
                % Only plot photons that exit at a certain angle.
                if (dirz < acceptance_angle)
                    continue;
                end
                
                % retrieve the stored path length in the medium.
                path_length = data(k,5);
                
                % retrieve the exit location of the photon.
                x = data(k,6);
                y = data(k,7);
                z = data(k,8);
                
                % add the distance from medium to the detector pixel.
                dist_to_pixel = sqrt(D^2 + (x_pixel - x)^2 + (y_pixel - y)^2);
                L = path_length + dist_to_pixel;
                
                CCDGrid(i,j) = (CCDGrid(i,j) + 1/(dist_to_pixel*weight*exp(-(1i)*L*2*pi/lambda)));
            end
        end
        y_pixel = start_y + CCDdy;
    end
    
    figure(speckleFigure);
    imagesc((abs(CCDGrid)))
    colormap hot;
    drawnow;
    if (createMovie)
        set(gcf,'Renderer', 'zbuffer');
        frame = getframe(gcf);
        aviobj = addframe(aviobj, frame);
    end
end

if (createMovie)
    aviobj = close(aviobj);
end




% -------------------------------------------------------------------------


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