

t_start = 1;
t_end   = 197;


aviobj = avifile ( 'speckle-modulation', 'compression', 'None', 'fps', 15);

figure;

for i=t_start:t_end
    
    
        speckleData = dlmread(['speckle-', num2str(i), '.txt']);
        % Normalize the data.
        speckleData = speckleData ./ max(max(max(speckleData)));
        imagesc(speckleData);
        %cmax = caxis;
        %caxis([cmax(2)/3 cmax(2)]);  % Scale the colormap of the image.
        %caxis auto;
        colormap hot;
        drawnow;
    
    
    set(gcf,'Renderer', 'zbuffer');
    frame = getframe(gcf);
    aviobj = addframe(aviobj, frame);
end