function varargout = GUI_ref_asc_desc(varargin)
% GUI_REF_ASC_DESC MATLAB code for GUI_ref_asc_desc.fig
%      GUI_REF_ASC_DESC, by itself, creates a new GUI_REF_ASC_DESC or raises the existing
%      singleton*.
%
%      H = GUI_REF_ASC_DESC returns the handle to a new GUI_REF_ASC_DESC or the handle to
%      the existing singleton*.
%
%      GUI_REF_ASC_DESC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_REF_ASC_DESC.M with the given input arguments.
%
%      GUI_REF_ASC_DESC('Property','Value',...) creates a new GUI_REF_ASC_DESC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_ref_asc_desc_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_ref_asc_desc_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_ref_asc_desc

% Last Modified by GUIDE v2.5 04-Dec-2020 16:00:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_ref_asc_desc_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_ref_asc_desc_OutputFcn, ...
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


% --- Executes just before GUI_ref_asc_desc is made visible.
function GUI_ref_asc_desc_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_ref_asc_desc (see VARARGIN)
global happy_var happy_desc happy_asc im seg_I centers_asc centers_desc happy_automatic ObjectPolarity_mode
happy_desc = 'no'; happy_asc = 'no'; happy_automatic = 'no'; centers_asc = []; centers_desc = [];
handles.output = hObject;
im = getappdata(0,'im')
seg_I = [];
ObjectPolarity_mode = getappdata(0,'ObjectPolarity_mode')
axes(handles.axes1), imshow(im)
set(handles.Happy_str, 'Visible', 'off');set(handles.Yes_button, 'Visible', 'off');set(handles.No_button, 'Visible', 'off');

% Choose default command line output for GUI_ref_asc_desc
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_ref_asc_desc wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_ref_asc_desc_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
global im centers_desc radii_desc centers_asc radii_asc xmin_desc xmax_desc ymin_desc ymax_desc xmin_asc xmax_asc ymin_asc ymax_asc happy_asc happy_desc happy_var happy_automatic
if ~isempty(centers_asc) && ~isempty(centers_desc)
    [xmin_asc xmax_asc ymin_asc ymax_asc xmin_desc xmax_desc ymin_desc ymax_desc] = find_frame_coordinates(im, centers_asc, centers_desc)
end
outputname = {'centers_desc','radii_desc','centers_asc','radii_asc','xmin_desc','xmax_desc','ymin_desc','ymax_desc','xmin_asc','xmax_asc','ymin_asc','ymax_asc'};
handles.output = [{centers_desc} radii_desc {centers_asc} radii_asc xmin_desc xmax_desc ymin_desc ymax_desc xmin_asc xmax_asc ymin_asc ymax_asc];
varargout{1} = handles.output;
for k=1:length(outputname)
    eval(['setappdata(0,''' outputname{k} ''',handles.output{k})'])
end
if strcmp(happy_automatic,'yes') || (strcmp(happy_asc,'yes') && strcmp(happy_desc,'yes'))
    uiresume
    close
end


% --- Executes on button press in Automatic_selection_button.
function Automatic_selection_button_Callback(hObject, eventdata, handles)
global im seg_I centers_desc radii_desc centers_asc radii_asc xmin_desc xmax_desc ymin_desc ymax_desc xmin_asc xmax_asc ymin_asc ymax_asc happy_asc happy_desc happy_var happy_automatic
happy_var = 'no';
if strcmp(happy_var,'no')
    if isempty(seg_I)
        [centers_desc radii_desc centers_asc radii_asc xmin_desc xmax_desc ymin_desc ymax_desc xmin_asc xmax_asc ymin_asc ymax_asc] = find_desc_asc_from_ref_scans(im);
        axes(handles.axes1),imshow(im), viscircles(centers_asc, radii_asc,'Color','r'), viscircles(centers_desc, radii_desc,'Color','b');
    else
        [centers_desc radii_desc centers_asc radii_asc xmin_desc xmax_desc ymin_desc ymax_desc xmin_asc xmax_asc ymin_asc ymax_asc] = find_desc_asc_from_ref_scans(seg_I);
        axes(handles.axes1),imshow(im), viscircles(centers_asc, radii_asc,'Color','r'), viscircles(centers_desc, radii_desc,'Color','b');
    end
    set(handles.Happy_str, 'Visible', 'on');set(handles.Yes_button, 'Visible', 'on');set(handles.No_button, 'Visible', 'on');
    waitforbuttonpress;
    pause(.4)
    if strcmp(happy_var,'yes')
        happy_automatic = 'yes'
        uiresume
    else
        pause(.4)
        centers_desc = [], centers_asc = [];
        happy_automatic = 'no'
        imshow(im)
        set(handles.Happy_str, 'Visible', 'off');set(handles.Yes_button, 'Visible', 'off');set(handles.No_button, 'Visible', 'off');
    end
end

%- Store variable in handle
guidata(hObject, handles);

% --- Executes on button press in Ascending_aorta_button.
function Ascending_aorta_button_Callback(hObject, eventdata, handles)
global im seg_I centers_desc radii_desc centers_asc radii_asc xmin_desc xmax_desc ymin_desc ymax_desc xmin_asc xmax_asc ymin_asc ymax_asc happy_asc happy_desc happy_var happy_automatic ObjectPolarity_mode
if isempty(seg_I)
    [centers, radii] = imfindcircles(im,[1 100],'ObjectPolarity',ObjectPolarity_mode);
else
    [centers, radii] = imfindcircles(seg_I,[1 100],'ObjectPolarity',ObjectPolarity_mode);
end
%- Delete too large radii
ind = find(radii>25);
centers(ind,:) = [];
radii(ind) = [];
%% Ascending aorta
happy_var = 'no', happy_asc = 'no'; 
while strcmp(happy_var,'no')
    centers_asc = [];
    if ~isempty(centers_desc)
        imshow(im), viscircles(centers, radii,'Color','r'), viscircles(centers_desc, radii_desc,'Color','b');
    else
        imshow(im), viscircles(centers, radii,'Color','r');
    end
    [x,y] = ginput(1);
    [x_asc y_asc radii_asc] = find_ref(centers,radii,x,y);
    if ~isempty(centers_desc)
        imshow(im), viscircles([x_asc y_asc], radii_asc,'Color','r'), viscircles(centers_desc, radii_desc,'Color','b');
    else
        imshow(im), viscircles([x_asc y_asc], radii_asc,'Color','r');
    end
    set(handles.Happy_str, 'Visible', 'on');set(handles.Yes_button, 'Visible', 'on');set(handles.No_button, 'Visible', 'on');
    waitforbuttonpress;
    pause(.4)
    if strcmp(happy_var,'yes')
        centers_asc = [x_asc y_asc];
        happy_asc = 'yes';
        set(handles.Happy_str, 'Visible', 'off');set(handles.Yes_button, 'Visible', 'off');set(handles.No_button, 'Visible', 'off');
        if strcmp(happy_desc,'yes')
            uiresume
        else
            break
        end
    else
        set(handles.Happy_str, 'Visible', 'off');set(handles.Yes_button, 'Visible', 'off');set(handles.No_button, 'Visible', 'off');
    end
end
happy_var = 'no';

%- Store variable in handle
guidata(hObject, handles);

% --- Executes on button press in Descending_aorta_button.
function Descending_aorta_button_Callback(hObject, eventdata, handles)
global im seg_I centers_desc radii_desc centers_asc radii_asc xmin_desc xmax_desc ymin_desc ymax_desc xmin_asc xmax_asc ymin_asc ymax_asc happy_asc happy_desc happy_var happy_automatic ObjectPolarity_mode
if isempty(seg_I)
    [centers, radii] = imfindcircles(im,[1 100],'ObjectPolarity',ObjectPolarity_mode);
else
    [centers, radii] = imfindcircles(seg_I,[1 100],'ObjectPolarity',ObjectPolarity_mode);
end
%- Delete too large radii
ind = find(radii>25)
centers(ind,:) = [];
radii(ind) = [];
%% Descending aorta
happy_var = 'no', happy_desc = 'no';
while strcmp(happy_var,'no')
    centers_desc = [];
    if ~isempty(centers_asc)
        imshow(im), viscircles(centers, radii,'Color','b'), viscircles(centers_asc, radii_asc,'Color','r');
    else
        imshow(im), viscircles(centers, radii,'Color','b');
    end
    [x,y] = ginput(1);
    [x_desc y_desc radii_desc] = find_ref(centers,radii,x,y);
    if ~isempty(centers_asc)
        imshow(im), viscircles([x_desc y_desc], radii_desc,'Color','b'), viscircles(centers_asc, radii_asc,'Color','r');
    else
        imshow(im), viscircles([x_desc y_desc], radii_desc,'Color','b');
    end
    set(handles.Happy_str, 'Visible', 'on');set(handles.Yes_button, 'Visible', 'on');set(handles.No_button, 'Visible', 'on');
    waitforbuttonpress;
    pause(.4)
    if strcmp(happy_var,'yes')
        happy_desc = 'yes';
        centers_desc = [x_desc y_desc];
        set(handles.Happy_str, 'Visible', 'off');set(handles.Yes_button, 'Visible', 'off');set(handles.No_button, 'Visible', 'off');
        if strcmp(happy_asc,'yes')
            uiresume
        else
            break
        end
    else
        set(handles.Happy_str, 'Visible', 'off');set(handles.Yes_button, 'Visible', 'off');set(handles.No_button, 'Visible', 'off');
    end
end
% happy_var = 'no';

%- Store variable in handle
guidata(hObject, handles);


% --- Executes on button press in Forego_Asc_Aorta_button.
function Forego_Asc_Aorta_button_Callback(hObject, eventdata, handles)
% hObject    handle to Forego_Asc_Aorta_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
centers_asc = []; radii_asc = [];
happy_asc = 'yes';


% --- Executes on button press in Yes_button.
function Yes_button_Callback(hObject, eventdata, handles)
global happy_var
happy_var = 'yes';

% --- Executes on button press in No_button.
function No_button_Callback(hObject, eventdata, handles)
global happy_var
happy_var = 'no';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- ADDITIONAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [centers_desc_ref radii_desc_ref centers_asc_ref radii_asc_ref xmin_desc xmax_desc ymin_desc ymax_desc xmin_asc xmax_asc ymin_asc ymax_asc] = find_desc_asc_from_ref_scans(X)
    global ObjectPolarity_mode
    I = X;
    [centers, radii] = imfindcircles(I,[1 100],'ObjectPolarity',ObjectPolarity_mode);
    %- Delete too large radii
    ind = find(radii>25);
    centers(ind,:) = [];
    radii(ind) = [];

    %% Descending aorta
    %- Define rectangle
    %- Y
    y_width = size(I,1)/5;
    ymid = size(I,1)*2/3;
%     y_width = size(I,1)/7.5;
%     ymid = size(I,1)/2;
    ymin_desc = ymid*1.25 - y_width; ymax_desc = ymid*1.25 + y_width;
    %- X
    x_width = size(I,2)/7.5;
    xmid = size(I,2)/2;
    xmin_desc = xmid - x_width; xmax_desc = xmid + x_width;
%     imshow(I),rectangle('Position',[xmin_desc ymin_desc xmax_desc-xmin_desc ymax_desc-ymin_desc])
    
    centers_desc_potential = centers;radii_desc_potential = radii;
    %- Select circles within the middle rectangle
    for k=size(centers_desc_potential,1):-1:1
        if centers_desc_potential(k,1)<xmin_desc || centers_desc_potential(k,1)>xmax_desc || centers_desc_potential(k,2)<ymin_desc || centers_desc_potential(k,2)>ymax_desc
            centers_desc_potential(k,:) = [];
            radii_desc_potential(k,:) = [];
        end
    end
%     imshow(I),rectangle('Position',[xmin_desc ymin_desc xmax_desc-xmin_desc ymax_desc-ymin_desc]),viscircles(centers_desc_potential,radii_desc_potential,'Color','b')
    
    %- Select descending aorta
    [radii_desc_ref imax] = max(radii_desc_potential);
    centers_desc_ref = centers_desc_potential(imax,:);
    
    %% Ascending aorta
    %- Define middle location
    %- Y
    y_width = size(I,1)/5;
    ymid = size(I,1)/4 
%     y_width = size(I,1)/7.5;
%     ymid = size(I,1)/2 - y_width/2;
    ymin_asc = ymid - y_width; ymax_asc = ymid + y_width;
    %- X
    x_width = size(I,2)/7.5;
    xmid = size(I,2)/2;
    xmin_asc = xmid - x_width; xmax_asc = xmid + x_width;
%     figure(),imshow(I),rectangle('Position',[xmin_asc ymin_asc xmax_asc-xmin_asc ymax_asc-ymin_asc])
    
    %- Select circles within the middle rectangle
    for k=size(centers,1):-1:1
        if centers(k,1)<xmin_asc || centers(k,1)>xmax_asc || centers(k,2)<ymin_asc || centers(k,2)>ymax_asc
            centers(k,:) = [];
            radii(k,:) = [];
        end
    end
%     imshow(I), viscircles(centers, radii,'Color','b'),rectangle('Position',[xmin_asc ymin_asc xmax_asc-xmin_asc ymax_asc-ymin_asc])
    
    %- Select ascending aorta
    [radii_asc_ref imax] = max(radii);
    centers_asc_ref = centers(imax,:);
    
function [x_ref y_ref rad_ref] = find_ref(centers,radii,x,y);
for k=1:size(centers)
   d(k) = sqrt((centers(k,1)-x)^2 + (centers(k,2)-y)^2);
end
[c imin] = min(d)
x_ref = centers(imin,1);
y_ref = centers(imin,2);
rad_ref = radii(imin);
    
function [xmin_asc xmax_asc ymin_asc ymax_asc xmin_desc xmax_desc ymin_desc ymax_desc] = find_frame_coordinates(im, centers_asc, centers_desc)
%% Ascending aorta
%- X
x_width = size(im,2)/7.5;
xmin_asc = centers_asc(1) - x_width/2; xmax_asc = centers_asc(1) + x_width/2;
%- Y
y_width = size(im,1)/7.5;
ymin_asc = centers_asc(2) - y_width/2; ymax_asc = centers_asc(2) + y_width/2;

%% Descending aorta
%- X
x_width = size(im,2)/7.5;
xmin_desc = centers_desc(1) - x_width/2; xmax_desc = centers_desc(1) + x_width/2;
%- Y
y_width = size(im,1)/7.5;
ymin_desc = centers_desc(2) - y_width/2; ymax_desc = centers_desc(2) + y_width/2;
