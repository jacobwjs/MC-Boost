%plots pressure, photon path, etc.

% Nx = 64*2+4;
% Ny = 64*2+4;
% Nz = 64*2+4;

set(0,'DefaultFigureRenderer','OpenGL');

% number of photons to plot.
num_photons = 10;

% open the file that contains the photon coordinates.
%fid = fopen('photon-paths.txt');

% open the pressure file that k-wave genrated.
%pressure = dlmread('pressure-at-25us.txt');

%pressure = dlmread('pressure-at-9mm-2x2grid-2p5MHz.txt');
%pressure = dlmread('pressure-at-18mm-2x2grid-2p5MHz.txt');

% reshape into a 64x64x64 matrix.
% NOTE: this assumes using default dimensions from the simulation.
%pressure = reshape(pressure, Nx, Ny, Nz);



g_size = 4; %[cm]

% figure;
% view(45,30); 
% 
% %colormap(color_map); 
% caxis([0 1]); 
% axis([0 g_size 0 g_size 0 g_size]);
% %grid on;  
% 
% % plot orientation line
% a = 2*ones(1,g_size);
% b = 1:1:g_size;
% c = a;
% plot3(a,b,c);
% hold on;


% plot the pressure
[x,y,z] = meshgrid(0:g_size/Nz:g_size-g_size/Nz);
%isosurface(x,y,z,pressure);

% view(3);
% camlight; 
% lighting gouraud;
% 
% set(gcf,'Color',[.5,.5,.5],'Renderer','zbuffer');
% set(gca,'Color','black','XColor','white', ...
% 	'YColor','white','ZColor','white');






% Create movie of pressure
 aviobj = avifile('pressure20.avi','compression','None');
for k=1:1:353
    pressure = reshape(sensor_data(:,k), Nx, Ny, Nz);
    p1 = patch(isosurface(pressure,.1), ...
   'FaceColor','none','EdgeColor','blue');
p2 = patch(isocaps(pressure, .1),'FaceColor','interp',...
	'EdgeColor','none');
    %caxis([0 1]); 
    %axis([0 g_size 0 g_size 0 g_size]);
colormap(gray(100))
camlight left; camlight; lighting gouraud
isonormals(pressure,p1)
view(290,20);
    
%     isosurface(x,y,z,pressure);
%     caxis([0 1]); 
%     axis([0 g_size 0 g_size 0 g_size]);
%     
%     view(3);
%     camlight; 
%     lighting gouraud;
% 
    set(gcf,'Color',[.5,.5,.5],'Renderer','zbuffer');
    set(gca,'Color','black','XColor','white', ...
	'YColor','white','ZColor','white');
%     view(290,20);
%     box on;
    
    F = getframe(gcf);
    aviobj = addframe(aviobj,F);
    clf;
end
aviobj = close(aviobj);



% % plot the photon paths
% hold all;   % hold all cycles through colors for plotting.  Allows to distinguish
%             % between each separate photon path.
% for i=1:num_photons
%      tline = fgetl(fid);
%      tline = str2num(tline);
%      a = reshape(tline, 3, size(tline,2)/3);
%      plot3(a(3,:), a(2,:), a(1,:));
% end
% 
%  fclose(fid);
% 
% % add the axes labels
% xlabel('x [cm]');
% ylabel('z [cm]');
% zlabel('y [cm]');
% 
% view(3);
% 
% set(gcf,'Color',[.5,.5,.5],'Renderer','zbuffer')
% set(gca,'Color','black','XColor','white', ...
% 	'YColor','white','ZColor','white');
% 
% view(290,20);
% box on;
% camlight; 
% lighting gouraud;