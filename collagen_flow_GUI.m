function varargout = collagen_flow_GUI(varargin)
% COLLAGEN_FLOW_GUI M-file for collagen_flow_GUI.fig
%      COLLAGEN_FLOW_GUI, by itself, creates a new COLLAGEN_FLOW_GUI or raises the existing
%      singleton*.
%
%      H = COLLAGEN_FLOW_GUI returns the handle to a new COLLAGEN_FLOW_GUI or the handle to
%      the existing singleton*.
%
%      COLLAGEN_FLOW_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COLLAGEN_FLOW_GUI.M with the given input arguments.
%
%      COLLAGEN_FLOW_GUI('Property','Value',...) creates a new COLLAGEN_FLOW_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before collagen_flow_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to collagen_flow_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help collagen_flow_GUI

% Last Modified by GUIDE v2.5 14-May-2014 12:22:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @collagen_flow_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @collagen_flow_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else

    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before collagen_flow_GUI is made visible.
function collagen_flow_GUI_OpeningFcn(hObject, eventdata, handles, varargin) 
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to collagen_flow_GUI (see VARARGIN)

    % Choose default command line output for collagen_flow_GUI
    handles.output = hObject;
    
    %make the color_code figure invisible
    set(handles.color_code,'Visible','off')
    
    %define an initial data_path
    handles.d_path='C:/'
    
    %make sure that the units are updated
    set(handles.edit_length_scale_unit,'String',{'m','mm','µm','nm'})
    handles.time_scale_unit = 60;
    handles.time_scale_str = 'hour';
            handles.length_scale_unit = 1;
            handles.length_scale_str = 'µm';
            
    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes collagen_flow_GUI wait for user response (see UIRESUME)
    % uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = collagen_flow_GUI_OutputFcn(hObject, eventdata, handles)  %#ok<*INUSL>
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;


    % --- Executes on slider movement.
    function slider1_Callback(hObject, eventdata, handles)
    % hObject    handle to slider1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'Value') returns position of slider
    %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    files=handles.files;
    d_path=handles.d_path;
        i=round(get(handles.slider1,'Value'));
        im=imread([d_path,'/',files(i).name]);     
        im2=imread([d_path,'/',files(i+1).name]);
        if length(size(im))>2
            handles.im=im(:,:,get(handles.listbox_color,'Value'));       
            handles.im2=im2(:,:,get(handles.listbox_color,'Value'));
        else
            handles.im=im;     
            handles.im2=im2;
        end


    handles=detect_edge(hObject, eventdata, handles)
    guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles) %#ok<*DEFNU,*INUSD>
    % hObject    handle to slider1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: slider controls usually have a light gray background.
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end


% --------------------------------------------------------------------
function load_im_Callback(hObject, eventdata, handles)
% hObject    handle to load_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function load_folder_Callback(hObject, eventdata, handles)
% hObject    handle to load_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    

%make the color_code figure invisible
set(handles.color_code,'Visible','off')

set(handles.slider1,'Value',1);
%select the previously chosen path
d_path=handles.d_path;
d_path=uigetdir(d_path);
    files=dir([d_path,'/*.png']);
    curr=round(get(handles.slider1,'Value'));
    i=curr;
    im=imread([d_path,'/',files(i).name]);
    im2=imread([d_path,'/',files(i+1).name]);
    if length(size(im))>2
        handles.im=im(:,:,get(handles.listbox_color,'Value'));      
        
        handles.im2=im2(:,:,get(handles.listbox_color,'Value'));
    else
        handles.im=im;     
        handles.im2=im2;
    end
    set(handles.slider1,'Max',length(files)-1);
    slider_step(1) = 1/(length(files)-1);
    slider_step(2) = 1/(length(files)-1);
    set(handles.slider1,'sliderstep',slider_step);
    

handles.num_images=length(files);
axes(handles.image);
imshow(handles.im);
handles.curr=curr;
handles.tiff=0;
handles.d_path=d_path;
handles.files=files;
handles=detect_edge(hObject, eventdata, handles);
guidata(hObject, handles);



% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=correlate_im(hObject, eventdata, handles);

%interpolate without time averging
flow_struct=save_actual_im(hObject, eventdata, handles);
interpolate_flow_nested_single(handles)

guidata(hObject, handles);


function s_source_Callback(hObject, eventdata, handles)
% hObject    handle to s_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of s_source as text
%        str2double(get(hObject,'String')) returns contents of s_source as a double


% --- Executes during object creation, after setting all properties.
function s_source_CreateFcn(hObject, eventdata, handles)
% hObject    handle to s_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function s_search_Callback(hObject, eventdata, handles)
% hObject    handle to s_search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of s_search as text
%        str2double(get(hObject,'String')) returns contents of s_search as a double


% --- Executes during object creation, after setting all properties.
function s_search_CreateFcn(hObject, eventdata, handles)
% hObject    handle to s_search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function g_dist_Callback(hObject, eventdata, handles)
% hObject    handle to g_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of g_dist as text
%        str2double(get(hObject,'String')) returns contents of g_dist as a double


% --- Executes during object creation, after setting all properties.
function g_dist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to g_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function c_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to c_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c_thresh as text
%        str2double(get(hObject,'String')) returns contents of c_thresh as a double


% --- Executes during object creation, after setting all properties.
function c_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%This is the case where we will run the full dataset at the end we will
%check if the interpolation is checked, and if so it will be run.

%turn off warning:
warning('all','off')
    files=handles.files;
    d_path=handles.d_path;
    num_images=handles.num_images;

for i=1:num_images-1
    set(handles.slider1,'value',i)
    if handles.tiff==0
        im=imread([d_path,'/',files(i).name]);     
        im2=imread([d_path,'/',files(i+1).name]);
    else
        im=imread(files,i); 
        im2=imread(files,i+1);
    end
    if length(size(im))>2
        handles.im=im(:,:,get(handles.listbox_color,'Value'));      
        handles.im2=im2(:,:,get(handles.listbox_color,'Value'));
    else
        handles.im=im;     
        handles.im2=im2;
    end
    
    handles.curr=i;
    guidata(hObject, handles);
    handles=detect_edge(hObject, eventdata, handles);
    handles=correlate_im(hObject, eventdata, handles);
    guidata(hObject, handles);
    flow_struct=save_actual_im(hObject, eventdata, handles);
end

%Now we check if the user want the interpolation, and if so we will run it.
if (get(handles.checkbox_,'Value')==1)
%     path=handles.d_path;
%     m2p=str2num(get(handles.edit_mue2pix,'String'));
%     fps=1/str2num(get(handles.edit_sec_per_frame,'String'));
%     xcorr_t=str2num(get(handles.c_thresh,'String'));
%     k_size=str2num(get(handles.edit_kernel_size,'String'));
%     k_sigma=str2num(get(handles.edit_kernel_sigma,'String'));
%     k_size_t=str2num(get(handles.edit_temp_kernel_size,'String'));
%     k_sigma_t=str2num(get(handles.edit_temp_sigma,'String'));
%     max_vel=str2num(get(handles.edit_max_flow_vel,'String'));
%     arrow_distance=str2num(get(handles.edit_arrow_distance,'String'));
    interpolate_flow_nested(handles)
end

%turn on warning:
warning('all','on')

function size_er_Callback(hObject, eventdata, handles)
% hObject    handle to size_er (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of size_er as text
%        str2double(get(hObject,'String')) returns contents of size_er as a double
handles=detect_edge(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function size_er_CreateFcn(hObject, eventdata, handles)
% hObject    handle to size_er (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function size_dil_Callback(hObject, eventdata, handles)
% hObject    handle to size_dil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of size_dil as text
%        str2double(get(hObject,'String')) returns contents of size_dil as a double

handles=detect_edge(hObject, eventdata, handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function size_dil_CreateFcn(hObject, eventdata, handles)
% hObject    handle to size_dil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

handles=detect_edge(hObject, eventdata, handles);
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function remove_edge_p_Callback(hObject, eventdata, handles)
% hObject    handle to remove_edge_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of remove_edge_p as text
%        str2double(get(hObject,'String')) returns contents of remove_edge_p as a double

handles=detect_edge(hObject, eventdata, handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function remove_edge_p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to remove_edge_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function handles=correlate_im(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gd=str2double(get(handles.g_dist,'String'));
gs=str2double(get(handles.s_source,'String'));
ss=str2double(get(handles.s_search,'String'));
c_th=str2double(get(handles.c_thresh,'String'));


curr=round(get(handles.slider1,'Value'));
axes(handles.image);
imshow(uint8(handles.im_e))
handles.curr;

im=handles.im;

im2=handles.im2;
BW2=handles.BW2;
se = strel('disk',gs);
BW2=imerode(BW2,se);
n=1;
for i=1:gd:size(im,1)
    for j=1:gd:size(im,2)
        if BW2(i,j)==1
            try 
                sub_roi=im(i-gs:i+gs,j-gs:j+gs);
            catch
                continue
            end
            try
                sub_area=im2(i-ss:i+ss,j-ss:j+ss);
            catch
                continue
            end
            try
                c = normxcorr2(sub_roi,sub_area);
            catch
                continue
            end
                [xroi,yroi]=size(sub_roi);
            [xarea,yarea]=size(sub_area);
            c=c(xroi:xarea,yroi:yarea);
            % figure, surf(c), shading flat
            % offset found by correlation
            [max_c, imax] = max(c(:));
            [ypeak, xpeak] = ind2sub(size(c),imax(1));
            [xc,yc]=size(c);
%             x_p(n)=(j);
%             y_p(n)=(i);
%             x_v(n)=xpeak-(xc+1)/2-1;
%             y_v(n)=ypeak-(yc+1)/2-1;
            x_vector(i,j)=xpeak-(xc+1)/2;
            y_vector(i,j)=ypeak-(yc+1)/2;
            corr_val(i,j)=max_c;
            n=n+1;
            
            if max_c>=c_th
                x_p(n)=(j);
                y_p(n)=(i);
                x_v(n)=xpeak-(xc+1)/2;
                y_v(n)=ypeak-(yc+1)/2;
                c_val(n)=max_c;
%                 axes(handles.color_code);
%                 surf(c)
%                 pause(.1)
                %max_c
            end
            
        end
    end
end
hold on
quiver(x_p,y_p,x_v,y_v,'r');

hold off

handles.x_pos=x_p;
handles.y_pos=y_p;
handles.x_vec=x_v;
handles.y_vec=y_v;
handles.c_val=c_val;

guidata(hObject, handles);
%now I create the grid


% --- Executes on slider movement.
function handles=detect_edge(hObject, eventdata, handles)
    % hObject    handle to slider1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'Value') returns position of slider
    %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    
    %%%
    % does something for edge detection, and does modify handle
    % but what the fuck does it do.
    
    %%%
    
    doEdgeDetect = get(handles.doEdgeDetect,'Value');
    
    curr=round(get(handles.slider1,'Value'));
    set(handles.text1,'String',['Image: ',num2str(curr)]);
    axes(handles.image);
    handles.curr=curr;
    im=handles.im;
    im=im;
    imshow(im);
    guidata(hObject, handles);

    if(doEdgeDetect == 1.0)
        axes(handles.color_code);
        level = graythresh(im);

        %size_g=55;%size of the gaussian kernel
        size_g=str2double(get(handles.size_er,'String'));
        sig=1;%sigma of the gaussian kernel
        d=10;%edge rechtalnges to compensate inhomogenous illumination
        %t=1;%should I invert?
        t=get(handles.invert_image,'Value');
        ones1=5;%first imclose, removes smal particle
        %ones2=50;%second imclose, removes particles at the end of procedure and clese the contour
        ones2=str2double(get(handles.size_dil,'String'));
        %min_size=50;%remove particles after threshold
        min_size=str2double(get(handles.remove_edge_p,'String'));

        axes(handles.color_code);
        %I want that center is bright, so inverse if necessary
        if t==1
            im=double(im);
            %im=uint16((im-2^16)*-1);
            im=uint8((im-2^8)*-1);


        end

        %Now I correct for assymtric illumination. I take the corners and
        %intermolate a map for normalisation

        [x,y]=size(im);
        I(1,1)=mean(reshape(im(1:d,1:d),1,d^2));
        I(2,2)=mean(reshape(im(x-d+1:x,y-d+1:y),1,d^2));
        I(2,1)=mean(reshape(im(1:d,y-d+1:y),1,d^2));
        I(1,2)=mean(reshape(im(x-d+1:x,1:d),1,d^2));

        %now I interpolare the normalisation matrix
        [X,Y] = meshgrid(1:x-1:x,1:y-1:y);
        Z = I;
        [XI,YI] = meshgrid(1:x,1:y);
        ZI = interp2(X,Y,Z,XI,YI);

        %now I normalize the image
        %make sure the image is 8bit
        im=uint8(im);
        im_n=im-(uint8(ZI')-min(min(ZI)));

        %# Create the gaussian filter with hsize = [5 5] and sigma = 2
        G = fspecial('gaussian',[size_g size_g],sig);
        %# Filter it
        im_blur = imfilter(im_n,G,'same');
        % now threshold
        graythresh(im_blur)
        im_th=im2bw(im_blur, graythresh(im_blur));

        %now I clean up
        bw2 = imfill(im_th,'holes');
        %bw3 = imopen(bw2, ones(ones1,ones1));
        bw3 = imopen(im_th, ones(ones1,ones1));
        bw4 = bwareaopen(bw3, min_size);
        bw4b=(imclose(bw4,ones(ones2,ones2)));
        bw5=(imfill(bw4b,'holes'));
        imshow(im_n)
        %finally remove all particles that don't take up half of the total white
        %pixels numbers
        bw6= bwareaopen(bw5, round(sum(sum(bw5))/2));

        bw6_perim = bwperim(bw6);
        overlay1 = imoverlay(im, bw6_perim, [1 0 0]);

        axes(handles.image);
        imshow(overlay1);
        handles.BW2=(bw6-1)*-1;
    else
        handles.BW2=im*0+1;
    end

    
    handles.im_e=im;

    guidata(hObject, handles);


function flow_struct=save_actual_im(hObject, eventdata, handles)

im=handles.im;
files=handles.files;
curr=round(get(handles.slider1,'value'));
%First I arrange the structure
flow_struct.d_path=handles.d_path;
flow_struct.curr=curr;
flow_struct.im=im;
flow_struct.x_pos=handles.x_pos;
flow_struct.y_pos=handles.y_pos;
flow_struct.x_vec=handles.x_vec;
flow_struct.y_vec=handles.y_vec;
flow_struct.c_val=handles.c_val;
flow_struct.file=files(curr).name;
flow_struct.edge_im=handles.BW2;

e_meth.gaussk=str2num(get(handles.size_er,'String'));
e_meth.connect=str2num(get(handles.size_dil,'String'));
e_meth.remove=str2num(get(handles.remove_edge_p,'String'));
e_meth.invert=get(handles.invert_image,'Value');


flow_struct.edge_meth=e_meth;

flow_struct.size_cut=str2num(get(handles.remove_edge_p,'String'));
flow_struct.se_size=str2num(get(handles.size_dil,'String'));
flow_struct.se2_size=str2num(get(handles.size_er,'String'));


flow_struct.gd=str2num(get(handles.g_dist,'String'));
flow_struct.gs=str2num(get(handles.s_source,'String'));
flow_struct.ss=str2num(get(handles.s_search,'String'));
flow_struct.c_th=str2num(get(handles.c_thresh,'String'));

[status,message,messageid] = mkdir(handles.d_path,'data');
[status,message,messageid] = mkdir(handles.d_path,'correlation_image');
%this is not done
fil=files(curr).name;
axes(handles.image);
save([handles.d_path,'/data/',fil(1:end-3),'mat'],'flow_struct')
%print([handles.d_path,'/image/',fil(1:end-4),'_results.png'],'-fhandles.image','-dpng','-r300')
%saveas(handles.image,[handles.d_path,'/image/',fil(1:end-4),'_results.png'])
newfigure = figure;
imshow(uint8(handles.im_e))
hold on
quiver(handles.x_pos,handles.y_pos,handles.x_vec,handles.y_vec,'r');
% print(gcf,[handles.d_path,'/image/',fil(1:end-4),'_results.png'],'-dpng','-r150')
saveas(gcf,[handles.d_path,'/correlation_image/',fil(1:end-4),'_results.png'])
% close the new figure
close(newfigure) 


% --- Executes on selection change in listbox_color.
function listbox_color_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_color contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_color


% --- Executes during object creation, after setting all properties.
function listbox_color_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kernel_ee_Callback(hObject, eventdata, handles)
% hObject    handle to kernel_ee (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kernel_ee as text
%        str2double(get(hObject,'String')) returns contents of kernel_ee as a double


% --- Executes during object creation, after setting all properties.
function kernel_ee_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kernel_ee (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kernel_size_Callback(hObject, eventdata, handles)
% hObject    handle to kernel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kernel_size as text
%        str2double(get(hObject,'String')) returns contents of kernel_size as a double


% --- Executes during object creation, after setting all properties.
function kernel_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kernel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles=interpolate_flow(hObject,handles)
%Here we will apply an gaussian interpolation. Each measured point will
%have a gaussian field associated with it. The field will have the decay
%and the size as defined in the GUI. In the end, the average of each values
%is calculated and the image field is displayed and stored in the handle.
%Later it will be also put into the save structure that is stored if
%desired.



function edit_mue2pix_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mue2pix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mue2pix as text
%        str2double(get(hObject,'String')) returns contents of edit_mue2pix as a double


% --- Executes during object creation, after setting all properties.
function edit_mue2pix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mue2pix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sec_per_frame_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sec_per_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sec_per_frame as text
%        str2double(get(hObject,'String')) returns contents of edit_sec_per_frame as a double


% --- Executes during object creation, after setting all properties.
function edit_sec_per_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sec_per_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_.
function checkbox__Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_



function edit_kernel_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_kernel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kernel_size as text
%        str2double(get(hObject,'String')) returns contents of edit_kernel_size as a double


% --- Executes during object creation, after setting all properties.
function edit_kernel_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kernel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_kernel_sigma_Callback(hObject, eventdata, handles)
% hObject    handle to edit_kernel_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kernel_sigma as text
%        str2double(get(hObject,'String')) returns contents of edit_kernel_sigma as a double


% --- Executes during object creation, after setting all properties.
function edit_kernel_sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kernel_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_temp_kernel_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_temp_kernel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_temp_kernel_size as text
%        str2double(get(hObject,'String')) returns contents of edit_temp_kernel_size as a double


% --- Executes during object creation, after setting all properties.
function edit_temp_kernel_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_temp_kernel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_temp_sigma_Callback(hObject, eventdata, handles)
% hObject    handle to edit_temp_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_temp_sigma as text
%        str2double(get(hObject,'String')) returns contents of edit_temp_sigma as a double


% --- Executes during object creation, after setting all properties.
function edit_temp_sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_temp_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_max_flow_vel_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max_flow_vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_max_flow_vel as text
%        str2double(get(hObject,'String')) returns contents of edit_max_flow_vel as a double


% --- Executes during object creation, after setting all properties.
function edit_max_flow_vel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max_flow_vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_arrow_distance_Callback(hObject, eventdata, handles)
% hObject    handle to edit_arrow_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_arrow_distance as text
%        str2double(get(hObject,'String')) returns contents of edit_arrow_distance as a double


% --- Executes during object creation, after setting all properties.
function edit_arrow_distance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_arrow_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%In this function we will start to do the interpolation
function interpolate_flow_nested(handles)
    %Here we will load the retrograde flow data, and then interpolate it
    %according to the spatial and temporal kernel size. The result will be
    %stored in a subfolder. There will be a filtering: If the direction of a
    %flow changes apruptly in time of space, it will be ignored.
    %path='E:/Science/data/brian_stemer/temp1'
    filter_versioninfo='advanced_display_and_filter_v1_5';
    path=handles.d_path;
    mue2pix_ratio=str2double(get(handles.edit_mue2pix,'String'));
    
    
    s = str2double(get(handles.edit_sec_per_frame,'String'));
    m = str2double(get(handles.min_per_frame,'String'));
    h = str2double(get(handles.hrs_per_frame,'String'));
    
    fps = 1/(s+60*m+3600*h);
    
    clear s m h;
    
    xcorr_thresh=str2double(get(handles.c_thresh,'String'));
    k_size=2*(ceil(0.5*str2double(get(handles.edit_kernel_size,'String'))/mue2pix_ratio ));
    k_sigma=2*(ceil(0.5*str2double(get(handles.edit_kernel_sigma,'String'))/mue2pix_ratio ));
    k_size_temp=str2double(get(handles.edit_temp_kernel_size,'String'));
    k_sigma_temp=str2double(get(handles.edit_temp_sigma,'String'));
    
    max_flow_vel=str2double(get(handles.edit_max_flow_vel,'String'));
    max_flow_vel = max_flow_vel*handles.length_scale_unit/handles.time_scale_unit;
    flow_field_arrow_distance=str2num(get(handles.edit_arrow_distance,'String'));
 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %Define kernel size. All kernel values have to be even numbers!!!!!!
    % if nargin < 4, xcorr_thresh=0.6; end                                                 %threshold value for cross-correlation
    % if nargin < 5, k_size=2*ceil(ceil(2*5/mue2pix_ratio)/2), end       %spatial gaussian filter kernel, the radius will be 5µm 
    % if nargin < 6, k_sigma=2*ceil(ceil(1/mue2pix_ratio)/2), end       %sigma of the spatial gaussian kernel should be 1µm
    % if nargin < 7, k_size_temp=4; end                                                    %time gaussian filter kernel, +- 2 frames
    % if nargin < 8, k_sigma_temp=2; end                                                %sigma of time gaussian
    % if nargin < 9, max_flow_vel=4; end                                                   %maximal flow velocity
    % if nargin < 10, flow_field_arrow_distance=25; end                               %distance between arrows in the interpolated images                    %maximal flow velocity
    % 
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %
    %Here we create the folder where we will later store the interpolation
    %images and matlab files
    mkdir([path,filesep,'interpolation_data']);


    %first we will read the file list:
    files=dir([path,filesep,'data',filesep,'*mat'])
    for i=1:length(files)
        load([path,filesep,'data',filesep,files(i).name]);
        flow_data(i)=flow_struct;
        clear list
        %We directly switch from pixel per image to the mue/min
        list(:,1)=flow_struct.x_pos;
        list(:,2)=flow_struct.y_pos;
        list(:,3)=flow_struct.x_vec*mue2pix_ratio*fps*60;
        list(:,4)=flow_struct.y_vec*mue2pix_ratio*fps*60;
        list(:,5)=flow_struct.c_val;
        list_cell{i}=list;
    end

    %get the max and min x and y values of the rf vector startpoints for the full image series (means
    %you get the values closest to the image boarder in all the time series)
    t=max(size(list_cell));
    for j=1:t
        list=list_cell{j};
        y(j)=max(list(:,2));
        x(j)=max(list(:,1));
        y_min(j)=min(list(:,2));
        x_min(j)=min(list(:,1));
    end

    %to prevent problems with rf values too close to the upper or left image boarder (so that the filter kernel would
    %reach outside the image and negative array indices are not defined, in the positive direction there is no problem, the code
    %just fills up the missing values with zeros),
    %check if the min x,y values of the retroflow points in list_cell are within kernel/2.
    %if not, we have to shift the startpoints of the flow arrows, and all the images in the series.
    %all the shift info is stored in the x_shift, and y_shift variable.
    if (min(x_min)<=k_size/2), x_shift=ceil(k_size/2-min(x_min)+1); else x_shift=0; end
    if (min(y_min)<=k_size/2), y_shift=ceil(k_size/2-min(y_min)+1); else y_shift=0; end

    %Save all parameters including x_shift and y_shift in a parameter file:
    save([path,filesep,'analysis_parameters.mat'],'filter_versioninfo','path','xcorr_thresh','k_size','k_sigma','k_size_temp','k_sigma_temp','mue2pix_ratio','x_shift','y_shift');

    %now shift the retroflow values for the full image series
    for j=1:t
        list=list_cell{j};
        list(:,2)=list(:,2)+y_shift;
        list(:,1)=list(:,1)+x_shift;
        list_cell{j}=list;    
    end

    %define the k_size_temp+1 (e.g 5) layer thick convolution array block that is moved through the full
    %image series in the analysis.
    x=max(x);
    y=max(y);
    im_x(1:y+k_size+1+y_shift,2:x+k_size+1+x_shift,1:k_size_temp+1)=0;
    im_devider=im_x;
    im_y=im_x;

    %generate the convolution kernel, gaussian kernel, circular, sigma of k_size/6
    [k_x,k_y]=meshgrid(-k_size/2:+k_size/2,-k_size/2:+k_size/2);
    kernel_2D=exp(-(k_x.^2+k_y.^2)/(2*(k_sigma)^2));
    for i=1:k_size_temp/2+1
        k_size_temp-i+2;
        kernel_3D(:,:,i)=kernel_2D*exp(-(k_size_temp/2+1-i)^2/(k_sigma_temp)^2);
        kernel_3D(:,:,k_size_temp-i+2)=kernel_3D(:,:,i);    
    end
    min_val=max(min(kernel_3D(:,:,k_size_temp/2+1)));
    kernel_3D(find(kernel_3D<min_val))=0;

    startframe=1;
    endframe=max(size(list_cell));
    %endframe=3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Now that we have the unfiltered retroflow data, we can start the filtering 
    %Do the 3D convolution
    for j=startframe:endframe;
        list=list_cell{j};
        j;
        tic
        for i=1:max(size(list));
            x_c=list(i,1);
            y_c=list(i,2);
            dx=list(i,3);
            dy=list(i,4);
    %        [th_d,r_d]=cart2pol(com_evolution(j,1)-x_c,com_evolution(j,2)-y_c);
            [th_f,r_f]=cart2pol(dx,dy);
    %        th_diff=abs(th_d-th_f)*360/(2*pi);
    %        if (th_diff>180), th_diff=abs(th_diff-360); end
            if(list(i,5)>xcorr_thresh) %&& th_diff<max_angle)
                x_c=list(i,1);
                y_c=list(i,2);
                %lb_y=y_c-(k_size/2);
                %ub_y=y_c+(k_size/2);
                %lb_x=x_c-(k_size/2);
                %ub_x=x_c+(k_size/2);
                %build the convolution matrix and the normalization matrix 
                im_x(y_c-(k_size/2):y_c+(k_size/2),x_c-(k_size/2):x_c+(k_size/2),:)=list(i,5)*list(i,3)*kernel_3D+im_x(y_c-(k_size/2):y_c+(k_size/2),x_c-(k_size/2):x_c+(k_size/2),:);
                im_y(y_c-(k_size/2):y_c+(k_size/2),x_c-(k_size/2):x_c+(k_size/2),:)=list(i,5)*list(i,4)*kernel_3D+im_y(y_c-(k_size/2):y_c+(k_size/2),x_c-(k_size/2):x_c+(k_size/2),:);
                im_devider(y_c-(k_size/2):y_c+(k_size/2),x_c-(k_size/2):x_c+(k_size/2),:)=list(i,5)*kernel_3D+im_devider(y_c-(k_size/2):y_c+(k_size/2),x_c-(k_size/2):x_c+(k_size/2),:);          
            end
        end
        im_x_act=single(im_x(:,:,1)./im_devider(:,:,1));
        im_y_act=single(im_y(:,:,1)./im_devider(:,:,1));
        %create image and save final values (starts at frame k_size_temp/2 + 1 because of time filtering)
        if j>startframe-1+k_size_temp/2
            i=j-k_size_temp/2;
            %create the image, plus get the filtered final results
            create_retro_flow_image_collagen(path,flow_data(j).im,flow_data(j).edge_im,i,im_x_act,im_y_act,x_shift,y_shift,mue2pix_ratio,fps,max_flow_vel,flow_field_arrow_distance, handles.files, handles.length_scale_str, handles.time_scale_str)
           % [final_res_filtered,min_pr,max_pr,max_r_flow]=advanced_create_il_v1_5(rf_path,i,im_x_act,im_y_act,final_res,x_shift,y_shift,mue2pix_ratio);
        end
        %delete first time layer of convolution array block and define next one
        im_x(:,:,1)=[];
        im_y(:,:,1)=[];
        im_devider(:,:,1)=[];
        im_x(:,:,k_size_temp+1)=0;
        im_y(:,:,k_size_temp+1)=0;
        im_devider(:,:,k_size_temp+1)=0;
        toc
        %create the last image save final values (has to be done differently because of time filtering)
        if j==max(size(list_cell));     %this means the last 2 frames run again but its easy this way
            for n=1:k_size_temp/2
                i=j-k_size_temp/2+n;
                im_x_act=single(im_x(:,:,n)./im_devider(:,:,n));
                im_y_act=single(im_y(:,:,n)./im_devider(:,:,n));
               % final_res=load([rf_path,'/final_results/',files_fin_res(i).name]);
                %create the image, plus get the filtered final results
                create_retro_flow_image_collagen(path,flow_data(j).im,flow_data(j).edge_im,i,im_x_act,im_y_act,x_shift,y_shift,mue2pix_ratio,fps,max_flow_vel,flow_field_arrow_distance, handles.files, handles.length_scale_str, handles.time_scale_str)
           
    %            [final_res_filtered,min_pr,max_pr,max_r_flow]=advanced_create_il_v1_5(rf_path,i,im_x_act,im_y_act,final_res,x_shift,y_shift,mue2pix_ratio);
            end
        end

    end


    


% --- Executes on button press in invert_image.
function invert_image_Callback(hObject, eventdata, handles)
% hObject    handle to invert_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of invert_image



function min_per_frame_Callback(hObject, eventdata, handles)
% hObject    handle to min_per_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_per_frame as text
%        str2double(get(hObject,'String')) returns contents of min_per_frame as a double


% --- Executes during object creation, after setting all properties.
function min_per_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_per_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hrs_per_frame_Callback(hObject, eventdata, handles)
% hObject    handle to hrs_per_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hrs_per_frame as text
%        str2double(get(hObject,'String')) returns contents of hrs_per_frame as a double


% --- Executes during object creation, after setting all properties.
function hrs_per_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hrs_per_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_length_scale_unit_Callback(hObject, eventdata, handles)   
    contents = cellstr(get(hObject,'String'));
    switch contents{get(hObject,'Value')}
        case 'm'
            handles.length_scale_unit = 1e6;
            handles.length_scale_str = 'm';
        case 'mm'
            handles.length_scale_unit = 1e3;
            handles.length_scale_str = 'mm';
        case 'µm'
            handles.length_scale_unit = 1;
            handles.length_scale_str = 'µm';
        case 'nm'
            handles.length_scale_unit = 1e-3;    
            handles.length_scale_str = 'nm';
        otherwise
            assert(false,'unknown lenghtscale')
    end
    guidata(hObject, handles);



function edit_time_scale_unit_Callback(hObject, eventdata, handles)
    contents = cellstr(get(hObject,'String'));
    switch contents{get(hObject,'Value')}
        case 'day'
            handles.time_scale_unit = 1440;
            handles.time_scale_str = 'day';
        case 'hour'
            handles.time_scale_unit = 60;
            handles.time_scale_str = 'hour';
        case 'minute'
            handles.time_scale_unit = 1;
            handles.time_scale_str = 'minute';
        case 'second'
            handles.time_scale_unit = 1/60;
            handles.time_scale_str = 'second';
        otherwise
            assert(false,'unknown timescale')
            
    end
    guidata(hObject, handles);


    
    
    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% --- Executes during object creation, after setting all properties.
function edit_time_scale_unit_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to edit_time_scale_unit (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    handles.time_scale_unit = 1440;
    handles.time_scale_str = 'day';
    guidata(hObject, handles);
    




% --- Executes during object creation, after setting all properties.
function edit_length_scale_unit_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to edit_length_scale_unit (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    handles.length_scale_unit = 1e6;
    handles.length_scale_str = 'm';
    guidata(hObject, handles);


% --- Executes on button press in doEdgeDetect.
function doEdgeDetect_Callback(hObject, eventdata, handles)
% hObject    handle to doEdgeDetect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of doEdgeDetect



%In this function we will start to do the interpolation
function interpolate_flow_nested_single(handles)
    %Here we will load the retrograde flow data, and then interpolate it
    %according to the spatial and temporal kernel size. The result will be
    %stored in a subfolder. There will be a filtering: If the direction of a
    %flow changes apruptly in time of space, it will be ignored.
    %path='E:/Science/data/brian_stemer/temp1'
    filter_versioninfo='advanced_display_and_filter_v1_5';
    path=handles.d_path;
    mue2pix_ratio=str2double(get(handles.edit_mue2pix,'String'));
    
    
    s = str2double(get(handles.edit_sec_per_frame,'String'));
    m = str2double(get(handles.min_per_frame,'String'));
    h = str2double(get(handles.hrs_per_frame,'String'));
    
    fps = 1/(s+60*m+3600*h);
    
    clear s m h;
    
    xcorr_thresh=str2double(get(handles.c_thresh,'String'));
    k_size=2*(ceil(0.5*str2double(get(handles.edit_kernel_size,'String'))/mue2pix_ratio ));
    k_sigma=2*(ceil(0.5*str2double(get(handles.edit_kernel_sigma,'String'))/mue2pix_ratio ));
    k_size_temp=0;
    k_sigma_temp=str2double(get(handles.edit_temp_sigma,'String'));
    
    max_flow_vel=str2double(get(handles.edit_max_flow_vel,'String'));
    max_flow_vel = max_flow_vel*handles.length_scale_unit/handles.time_scale_unit;
    flow_field_arrow_distance=str2num(get(handles.edit_arrow_distance,'String'));
 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %Define kernel size. All kernel values have to be even numbers!!!!!!
    % if nargin < 4, xcorr_thresh=0.6; end                                                 %threshold value for cross-correlation
    % if nargin < 5, k_size=2*ceil(ceil(2*5/mue2pix_ratio)/2), end       %spatial gaussian filter kernel, the radius will be 5µm 
    % if nargin < 6, k_sigma=2*ceil(ceil(1/mue2pix_ratio)/2), end       %sigma of the spatial gaussian kernel should be 1µm
    % if nargin < 7, k_size_temp=4; end                                                    %time gaussian filter kernel, +- 2 frames
    % if nargin < 8, k_sigma_temp=2; end                                                %sigma of time gaussian
    % if nargin < 9, max_flow_vel=4; end                                                   %maximal flow velocity
    % if nargin < 10, flow_field_arrow_distance=25; end                               %distance between arrows in the interpolated images                    %maximal flow velocity
    % 
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %
    %Here we create the folder where we will later store the interpolation
    %images and matlab files
    mkdir([path,filesep,'interpolation_data']);


    %first we will read the file list:
    %files=dir([path,filesep,'data',filesep,'*mat'])
    %for i=1:length(files)
    files=handles.files
    i=handles.curr
        load([path,filesep,'data',filesep,files(i).name(1:end-4),'.mat']);
        flow_data(1)=flow_struct;
        clear list
        %We directly switch from pixel per image to the mue/min
        list(:,1)=flow_struct.x_pos;
        list(:,2)=flow_struct.y_pos;
        list(:,3)=flow_struct.x_vec*mue2pix_ratio*fps*60;
        list(:,4)=flow_struct.y_vec*mue2pix_ratio*fps*60;
        list(:,5)=flow_struct.c_val;
        list_cell{1}=list;
   % end

    %get the max and min x and y values of the rf vector startpoints for the full image series (means
    %you get the values closest to the image boarder in all the time series)
    t=max(size(list_cell));
    for j=1:t
        list=list_cell{j};
        y(j)=max(list(:,2));
        x(j)=max(list(:,1));
        y_min(j)=min(list(:,2));
        x_min(j)=min(list(:,1));
    end

    %to prevent problems with rf values too close to the upper or left image boarder (so that the filter kernel would
    %reach outside the image and negative array indices are not defined, in the positive direction there is no problem, the code
    %just fills up the missing values with zeros),
    %check if the min x,y values of the retroflow points in list_cell are within kernel/2.
    %if not, we have to shift the startpoints of the flow arrows, and all the images in the series.
    %all the shift info is stored in the x_shift, and y_shift variable.
    if (min(x_min)<=k_size/2), x_shift=ceil(k_size/2-min(x_min)+1); else x_shift=0; end
    if (min(y_min)<=k_size/2), y_shift=ceil(k_size/2-min(y_min)+1); else y_shift=0; end

    %Save all parameters including x_shift and y_shift in a parameter file:
    save([path,filesep,'analysis_parameters.mat'],'filter_versioninfo','path','xcorr_thresh','k_size','k_sigma','k_size_temp','k_sigma_temp','mue2pix_ratio','x_shift','y_shift');

    %now shift the retroflow values for the full image series
    for j=1:t
        list=list_cell{j};
        list(:,2)=list(:,2)+y_shift;
        list(:,1)=list(:,1)+x_shift;
        list_cell{j}=list;    
    end

    %define the k_size_temp+1 (e.g 5) layer thick convolution array block that is moved through the full
    %image series in the analysis.
    x=max(x);
    y=max(y);
    im_x(1:y+k_size+1+y_shift,2:x+k_size+1+x_shift,1:k_size_temp+1)=0;
    im_devider=im_x;
    im_y=im_x;

    %generate the convolution kernel, gaussian kernel, circular, sigma of k_size/6
    [k_x,k_y]=meshgrid(-k_size/2:+k_size/2,-k_size/2:+k_size/2);
    kernel_2D=exp(-(k_x.^2+k_y.^2)/(2*(k_sigma)^2));
    for i=1:k_size_temp/2+1
        k_size_temp-i+2;
        kernel_3D(:,:,i)=kernel_2D*exp(-(k_size_temp/2+1-i)^2/(k_sigma_temp)^2);
        kernel_3D(:,:,k_size_temp-i+2)=kernel_3D(:,:,i);    
    end
    min_val=max(min(kernel_3D(:,:,k_size_temp/2+1)));
    kernel_3D(find(kernel_3D<min_val))=0;

    startframe=1;
    endframe=max(size(list_cell));
    %endframe=3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Now that we have the unfiltered retroflow data, we can start the filtering 
    %Do the 3D convolution
    for j=startframe:endframe;
        list=list_cell{j};
        j;
        tic
        for i=1:max(size(list));
            x_c=list(i,1);
            y_c=list(i,2);
            dx=list(i,3);
            dy=list(i,4);
    %        [th_d,r_d]=cart2pol(com_evolution(j,1)-x_c,com_evolution(j,2)-y_c);
            [th_f,r_f]=cart2pol(dx,dy);
    %        th_diff=abs(th_d-th_f)*360/(2*pi);
    %        if (th_diff>180), th_diff=abs(th_diff-360); end
            if(list(i,5)>xcorr_thresh) %&& th_diff<max_angle)
                x_c=list(i,1);
                y_c=list(i,2);
                %lb_y=y_c-(k_size/2);
                %ub_y=y_c+(k_size/2);
                %lb_x=x_c-(k_size/2);
                %ub_x=x_c+(k_size/2);
                %build the convolution matrix and the normalization matrix 
                im_x(y_c-(k_size/2):y_c+(k_size/2),x_c-(k_size/2):x_c+(k_size/2),:)=list(i,5)*list(i,3)*kernel_3D+im_x(y_c-(k_size/2):y_c+(k_size/2),x_c-(k_size/2):x_c+(k_size/2),:);
                im_y(y_c-(k_size/2):y_c+(k_size/2),x_c-(k_size/2):x_c+(k_size/2),:)=list(i,5)*list(i,4)*kernel_3D+im_y(y_c-(k_size/2):y_c+(k_size/2),x_c-(k_size/2):x_c+(k_size/2),:);
                im_devider(y_c-(k_size/2):y_c+(k_size/2),x_c-(k_size/2):x_c+(k_size/2),:)=list(i,5)*kernel_3D+im_devider(y_c-(k_size/2):y_c+(k_size/2),x_c-(k_size/2):x_c+(k_size/2),:);          
            end
        end
        im_x_act=single(im_x(:,:,1)./im_devider(:,:,1));
        im_y_act=single(im_y(:,:,1)./im_devider(:,:,1));
        %create image and save final values (starts at frame k_size_temp/2 + 1 because of time filtering)
        if j>startframe-1+k_size_temp/2
            i=j-k_size_temp/2;
            %create the image, plus get the filtered final results
            create_retro_flow_image_collagen(path,flow_data(j).im,flow_data(j).edge_im,handles.curr,im_x_act,im_y_act,x_shift,y_shift,mue2pix_ratio,fps,max_flow_vel,flow_field_arrow_distance, handles.files, handles.length_scale_str, handles.time_scale_str)
           % [final_res_filtered,min_pr,max_pr,max_r_flow]=advanced_create_il_v1_5(rf_path,i,im_x_act,im_y_act,final_res,x_shift,y_shift,mue2pix_ratio);
        end
        %delete first time layer of convolution array block and define next one
        im_x(:,:,1)=[];
        im_y(:,:,1)=[];
        im_devider(:,:,1)=[];
        im_x(:,:,k_size_temp+1)=0;
        im_y(:,:,k_size_temp+1)=0;
        im_devider(:,:,k_size_temp+1)=0;
        toc
        %create the last image save final values (has to be done differently because of time filtering)
        if j==max(size(list_cell));     %this means the last 2 frames run again but its easy this way
            for n=1:k_size_temp/2
                i=j-k_size_temp/2+n;
                im_x_act=single(im_x(:,:,n)./im_devider(:,:,n));
                im_y_act=single(im_y(:,:,n)./im_devider(:,:,n));
               % final_res=load([rf_path,'/final_results/',files_fin_res(i).name]);
                %create the image, plus get the filtered final results
                create_retro_flow_image_collagen(path,flow_data(j).im,flow_data(j).edge_im,handles.curr,im_x_act,im_y_act,x_shift,y_shift,mue2pix_ratio,fps,max_flow_vel,flow_field_arrow_distance, handles.files, handles.length_scale_str, handles.time_scale_str)
           
    %            [final_res_filtered,min_pr,max_pr,max_r_flow]=advanced_create_il_v1_5(rf_path,i,im_x_act,im_y_act,final_res,x_shift,y_shift,mue2pix_ratio);
            end
        end

    end
    %now I will plot the result in the right image on the GUI
    out_ca=[path,filesep,'color_code_and_arrows',filesep, files(handles.curr).name(1:end-4),'_','colorcode_arrows.png']
    im_cc=imread(out_ca);
    axes(handles.color_code)
    imshow(im_cc)
    axes(handles.image)

