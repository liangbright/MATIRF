function ShowDetectionResult(NameStr, DetectionData, ImageArray, Imin, Imax, Scale, ImageSaveInfo, PauseInterval)

%  Create and then hide the GUI as it is being constructed.

%  Construct the components of panel.
Panel_H = figure('name', NameStr, 'Visible','off','Position',[360,500,450,30],...
                 'CloseRequestFcn',{@CloseAll_Callback} );

Button_Play = uicontrol(Panel_H,...
                        'Style','pushbutton','String','Play',...
                        'Position',[1,1,70,25],...
                        'Callback',{@Play_Callback});
                         
Button_Stop = uicontrol(Panel_H,...
                       'Style','pushbutton','String','Stop',...
                       'Position',[80,1,70,25],...
                       'Callback',{@Stop_Callback});
                         
Button_Play_Stop_State=0; % 1 is play
                          % 0 is stop  

Button_Forward = uicontrol(Panel_H,...
                           'Style','pushbutton','String','Forward',...
                           'Position',[160,1,70,25],...
                           'Callback',{@Forward_Callback});
      
Button_Backward = uicontrol(Panel_H,...
                            'Style','pushbutton','String','Backword',...
                            'Position',[240,1,70,25],...
                            'Callback',{@Backward_Callback}); 

Button_SaveImage = uicontrol(Panel_H,...
                             'Style','pushbutton','String','SaveImage',...
                             'Position',[320,1,70,25],...
                             'Callback',{@SaveImage_Callback});                        

align([Button_Play, Button_Stop, Button_Forward, Button_Backward, Button_SaveImage],'fixed', 'Bottom');

% Change units to normalized so components resize automatically.
set([Panel_H, Button_Play, Button_Stop, Button_Forward, Button_Backward, Button_SaveImage],'Units','normalized');

%  Construct the image display area
Image_DispH =  figure('name', 'Image Frame', 'Visible','off');

% Initialize the GUI.
% Move the GUI to the center of the screen.
movegui(Panel_H,'north')
% Make the GUI visible.
set(Panel_H,'Visible','on');

%get image related parameters
[Ly Lx ImageNum]=size(ImageArray);

Pos=[500 500 Scale*Lx Scale*Ly];
set(Image_DispH, 'Position', Pos);
%movegui(Image_DispH,'south');

% Initialize the data
CurImageIdx=1;

Draw(1)

%  Callbacks for simple_gui. These callbacks automatically
%  have access to component handles and initialized data 
%  because they are nested at a lower level.

% Push button callbacks.
function Play_Callback(source,eventdata) 

if Button_Play_Stop_State ==0 
    Button_Play_Stop_State=1;
    
    set(Button_Play, 'enable','off');
    
    set(Button_Stop, 'enable','on');
    set(Button_Backward,'enable','off');
    set(Button_Forward,'enable','off');
    set(Button_SaveImage,'enable','off');
    
    if CurImageIdx >= ImageNum
        CurImageIdx =1;
    end        
end

while ishandle(Panel_H) && Button_Play_Stop_State == 1
         
    %{
    Draw(CurImageIdx);
    pause(PauseInterval);
    
    CurImageIdx=CurImageIdx+1;
    if CurImageIdx > ImageNum    
        CurImageIdx=1;
    end
    %}
       
    if CurImageIdx < ImageNum
       
        CurImageIdx=CurImageIdx+1;
        
        %display the image and the tracks in the axe.    
        Draw(CurImageIdx);
        
        pause(PauseInterval);
    else
        
        CurImageIdx=ImageNum;
        
        Button_Play_Stop_State=0;
        set(Button_Play, 'enable','on'); 
        set(Button_Stop, 'enable','off');  
        set(Button_Backward,'enable','on');
        set(Button_Forward,'enable','on');
        set(Button_SaveImage,'enable','on');
    end
    
end
    
end
 
 
%----------------------------------------
function Stop_Callback(source,eventdata) 

if Button_Play_Stop_State ==1

    Button_Play_Stop_State=0;

    set(Button_Stop, 'enable','off'); 
    
    set(Button_Play, 'enable','on');  
    set(Button_Backward,'enable','on');
    set(Button_Forward,'enable','on');
    set(Button_SaveImage,'enable','on');
end

end

%----------------------------------------
function Forward_Callback(source,eventdata) 

if CurImageIdx < ImageNum
       
    CurImageIdx=CurImageIdx+1;

else
    
    CurImageIdx=1;
    
end

%display the image and the tracks in the axe.
Draw(CurImageIdx);      

end
 
function Backward_Callback(source,eventdata) 

if CurImageIdx > 1       
    CurImageIdx=CurImageIdx-1;               
else
    CurImageIdx=ImageNum;
end
Draw(CurImageIdx);

end  


function Draw(t)
   
%display the image and the tracks in the axe.

I=ImageArray(:,:,t);

if ~ishandle(Image_DispH)
    figure(Image_DispH);
    Pos=[1 1 Scale*Lx Scale*Ly];
    set(Image_DispH, 'Position', Pos);    
end

%Pos = get(Image_DispH,'Position');
%Pos=[Pos(1) Pos(2) Scale*Lx Scale*Ly];
%set(Image_DispH, 'Position', Pos);

figure(Image_DispH); clf(Image_DispH); imagesc([],[], I, [Imin Imax]); 
colormap(gray(256)); axis off;
daspect([1 1 1])

ImageAxes=get(Image_DispH, 'CurrentAxes');
 
if DetectionData(t).ObsvNum>0

    theta=0:pi/30:2*pi;
    temp_x=cos(theta);
    temp_y=sin(theta);
    
    x=DetectionData(t).ObsvData(1,:);
    y=DetectionData(t).ObsvData(2,:);
    Sigma=DetectionData(t).ObsvData(4,:);  

    hold on
    %plot(ImageAxes, x, y,'r.');  
    
    ParticleNum=length(x);
    for k=1:ParticleNum     
        hold on    
        plot(ImageAxes, x(k)+temp_x*Sigma(k), y(k)+temp_y*Sigma(k),'g-');      
    end
end

drawnow;                      

impixelinfo(Image_DispH);

title(ImageAxes,['Image: ', num2str(t)]);

end

function SaveImage_Callback(source,eventdata)

set(Button_Play, 'enable','off'); 
set(Button_Stop, 'enable','off');  
set(Button_Backward,'enable','off');
set(Button_Forward,'enable','off');
set(Button_SaveImage,'enable','off');

for t=1:ImageNum
    Draw(t);        
    ImageAxes=get(Image_DispH, 'CurrentAxes');
    frame=getframe(ImageAxes);    
    tempI=frame2im(frame);    
    filename=[ImageSaveInfo.ImageDataDir FBSlashString ImageSaveInfo.FilenameBase index2str(t, 4) '.tif'];    
    imwrite(tempI, filename,'tif')    
    pause(0.1); 
end

set(Button_Play, 'enable','on'); 
set(Button_Stop, 'enable','on');  
set(Button_Backward,'enable','on');
set(Button_Forward,'enable','on');
set(Button_SaveImage,'enable','on');

end


function CloseAll_Callback(source,eventdata)

try    
    if ishandle(Image_DispH)    
        close(Image_DispH);        
    end
catch
    a=1;
end

delete(Panel_H)

end

% Note: keep the 'end' to make all the variables accessble to all the sub-functions
end
