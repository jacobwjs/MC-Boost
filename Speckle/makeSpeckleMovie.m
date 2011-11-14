

t_start = 1;
t_end   = 8;


%aviobj = avifile ( 'speckle-modulation', 'compression', 'None', 'fps', 2);
%figure;


F=moviein(t_end);
for i=t_start:t_end
    
    
    speckleData = dlmread(['speckle-', num2str(i), '.txt']);
    % Normalize the data.
    speckleData = speckleData ./ max(max(max(speckleData)));
    imagesc(speckleData);
    %cmax = caxis;
    %caxis([cmax(2)/3 cmax(2)]);  % Scale the colormap of the image.
    caxis auto;
    colormap hot;
    drawnow;
    
    F(:,i) = getframe;
    
    %{
    set(gcf,'Renderer', 'zbuffer');
    frame = getframe(gcf);
    aviobj = addframe(aviobj, frame);
    %}
    
end
% use 1st frame to get dimensions
[h, w, p] = size(F(1).cdata);
hf = figure; 
% resize figure based on frame's w x h, and place at (150, 150)
set(hf, 'Position', [150 150 w h]);
axis off
% Place frames at bottom left
repeat = 1;
fps = 15;
movie(hf,F,repeat,fps,[0 0 0 0]);
%aviobj = close(aviobj);