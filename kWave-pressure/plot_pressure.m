%plots pressure, photon path, etc.


set(0,'DefaultFigureRenderer','OpenGL');

% number of photons to plot.
num_photons = 10;

% open the file that contains the photon coordinates.
fid = fopen('../Default/photon-paths.txt');

% open the pressure file that k-wave genrated.
pressure = dlmread('pressure-at-25us.txt');

% reshape into a 64x64x64 matrix.
% NOTE: this assumes using default dimensions from the simulation.
pressure = reshape(pressure, 64, 64, 64);





figure;
view(45,30); 

%colormap(color_map); 
caxis([0 1]); 
axis([0 10 0 10 0 10]);
%grid on;  

% plot orientation line
a = 5*ones(1,10);
b = 1:1:10;
c = a;
plot3(a,b,c);
hold on;


% plot the pressure
[x,y,z] = meshgrid(0:10/64:10-10/64);
isosurface(x,y,z,pressure);

% plot the photon paths
hold all;   % hold all cycles through colors for plotting.  Allows to distinguish
            % between each separate photon path.
for i=1:num_photons
     tline = fgetl(fid);
     tline = str2num(tline);
     a = reshape(tline, 3, size(tline,2)/3);
     plot3(a(3,:), a(2,:), a(1,:));
end

 fclose(fid);

% add the axes labels
xlabel('x [cm]');
ylabel('z [cm]');
zlabel('y [cm]');

set(gcf,'Color',[.5,.5,.5],'Renderer','zbuffer')
set(gca,'Color','black','XColor','white', ...
	'YColor','white','ZColor','white');
box on;
camlight; lighting phong;