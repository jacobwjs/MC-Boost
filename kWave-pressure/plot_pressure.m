%plots pressure, photon path, etc.

p = dlmread('../Default/photon-paths.txt');
p = reshape(p, 3, size(p,2)/3);
[x,y,z] = meshgrid(0:10/64:10-10/64);

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
isosurface(x,y,z,pressure);

% plot the photon
plot3(p(3,:), p(2,:), p(1,:), '-m');

% add the axes labels
xlabel('x [cm]');
ylabel('z [cm]');
zlabel('y [cm]');

set(gcf,'Color',[.5,.5,.5],'Renderer','zbuffer')
set(gca,'Color','black','XColor','white', ...
	'YColor','white','ZColor','white');
box on;