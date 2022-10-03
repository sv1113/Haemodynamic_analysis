function varargout = GUI_ref_asc_desc_manual_axial_arch(varargin)
% GUI_REF_ASC_DESC_MANUAL_axial_arch MATLAB code for GUI_ref_asc_desc_manual_axial_arch.fig
%      _axial_arch, by itself, creates a new GUI_REF_ASC_DESC_MANUAL or raises the existing
%      singleton*.
%
%      H = GUI_REF_ASC_DESC_MANUAL_axial_arch returns the handle to a new GUI_REF_ASC_DESC_MANUAL or the handle to
%      the existing singleton*.
%
%      GUI_REF_ASC_DESC_MANUAL_axial_arch('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_REF_ASC_DESC_MANUAL.M with the given input arguments.
%
%      GUI_REF_ASC_DESC_MANUAL_axial_arch('Property','Value',...) creates a new GUI_REF_ASC_DESC_MANUAL_axial_arch or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_ref_asc_desc_manual_axial_arch_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_ref_asc_desc_manual_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help
% GUI_ref_asc_desc_manual

% Last Modified by GUIDE v2.5 28-Sep-2022 16:55:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_ref_asc_desc_manual_axial_arch_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_ref_asc_desc_manual_axial_arch_OutputFcn, ...
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
function GUI_ref_asc_desc_manual_axial_arch_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_ref_asc_desc_manual (see VARARGIN)
global happy_var happy_desc happy_asc im seg_I centers_asc centers_desc happy_automatic ObjectPolarity_mode flag_quit ...
    centers_desc_ref radii_desc_ref centers_asc_ref radii_asc_ref
happy_desc = 'no'; happy_automatic = 'no'; centers_asc = []; centers_desc = []; flag_quit = 0;
handles.output = hObject;
happy_asc = getappdata(0,'happy_asc')
im = getappdata(0,'im')
seg_I = [];
ObjectPolarity_mode = getappdata(0,'ObjectPolarity_mode');
ObjectPolarity_mode = getappdata(0,'ObjectPolarity_mode');
centers_desc_ref = getappdata(0,'centers_desc_ref');
radii_desc_ref = getappdata(0,'radii_desc_ref');
centers_asc_ref = getappdata(0,'centers_asc_ref');
radii_asc_ref = getappdata(0,'radii_asc_ref');

axes(handles.axes1), imshow(im)
set(handles.Happy_str, 'Visible', 'off');set(handles.Yes_button, 'Visible', 'off');set(handles.No_button, 'Visible', 'off');

% Choose default command line output for GUI_ref_asc_desc_manual
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_ref_asc_desc_manual wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_ref_asc_desc_manual_axial_arch_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
global im centers_desc radii_desc centers_asc radii_asc xmin_desc xmax_desc ymin_desc ymax_desc xmin_asc xmax_asc ymin_asc ymax_asc happy_asc happy_desc happy_var happy_automatic flag_quit
if ~isempty(centers_asc) && ~isempty(centers_desc)
    [xmin_asc xmax_asc ymin_asc ymax_asc xmin_desc xmax_desc ymin_desc ymax_desc] = find_frame_coordinates(im, centers_asc, centers_desc)
end

if strcmp(happy_asc,'yes') && radii_asc==0
    centers_asc = [0 0];
    radii_asc = 0;
end

outputname = {'centers_desc','radii_desc','centers_asc','radii_asc'};
handles.output = [{centers_desc} radii_desc {centers_asc} radii_asc];
varargout{1} = handles.output;

if flag_quit==1
    setappdata(0,'flag_quit',1)
    uiresume
elseif flag_quit==2
    setappdata(0,'flag_quit',2)
    uiresume
else
    setappdata(0,'flag_quit',0)
end

for k=1:length(outputname)
    eval(['setappdata(0,''' outputname{k} ''',handles.output{k})'])
end
if strcmp(happy_automatic,'yes') || (strcmp(happy_asc,'yes') && strcmp(happy_desc,'yes'))
    uiresume
    close
end


% --- Executes on button press in Automatic_selection_button.
function Automatic_selection_button_Callback(hObject, eventdata, handles)
global im centers_desc radii_desc centers_asc radii_asc  centers_desc_ref radii_desc_ref centers_asc_ref radii_asc_ref  ...
    happy_desc happy_var happy_automatic flag_quit
happy_var = 'no';
if strcmp(happy_var,'no')
    [centers_desc radii_desc centers_asc radii_asc] = find_desc_asc_from_ref_scans(im,centers_desc_ref,radii_desc_ref,centers_asc_ref,radii_asc_ref);
    flag_asc = getappdata(0,'flag_asc');
    axes(handles.axes1),imshow(im), viscircles(centers_desc, radii_desc,'Color','b');
    if flag_asc == 1
        viscircles(centers_asc, radii_asc,'Color','r')
    else
        centers_asc = [0 0]; radii_asc = 0;
    end
    set(handles.Happy_str, 'Visible', 'on');set(handles.Yes_button, 'Visible', 'on');set(handles.No_button, 'Visible', 'on');
    waitforbuttonpress;
    pause(.4)
    if strcmp(happy_var,'yes')
        happy_automatic = 'yes'
        uiresume
    else
%         pause(.4)
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
global im seg_I centers_desc radii_desc centers_asc radii_asc xmin_desc xmax_desc ymin_desc ymax_desc xmin_asc xmax_asc ymin_asc ymax_asc happy_asc happy_desc happy_var happy_automatic ObjectPolarity_mode flag_quit;
%% Ascending aorta
happy_var = 'no', happy_asc = 'no'; xo = []; yo = [];
while strcmp(happy_var,'no')
    centers_asc = [];
    imshow(im)
    if ~isempty(centers_desc)
        viscircles(centers_desc, radii_desc,'Color','b');
    end
    set(gcf,'windowbuttondownfcn',{@starttrack});
    waitfor(0,'userdata')
    asc_boundaries=get(0,'userdata');
    x_asc_boundaries = asc_boundaries(:,1); 
    y_asc_boundaries = handles.axes1.YLim(2)-asc_boundaries(:,2); 
    [xo yo radii_asc] = circle_fit(x_asc_boundaries,y_asc_boundaries);
    viscircles([xo yo], radii_asc,'Color','r')
        
    set(handles.Happy_str, 'Visible', 'on');set(handles.Yes_button, 'Visible', 'on');set(handles.No_button, 'Visible', 'on');
    waitforbuttonpress;
    pause(.4)
    if strcmp(happy_var,'yes')
        centers_asc = [xo yo];
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
global im seg_I centers_desc radii_desc centers_asc radii_asc xmin_desc xmax_desc ymin_desc ymax_desc xmin_asc xmax_asc ymin_asc ymax_asc happy_asc happy_desc happy_var happy_automatic ObjectPolarity_mode flag_quit
%% Descending aorta
happy_var = 'no', happy_desc = 'no'; xo = []; yo = [];
flag_asc = getappdata(0,'flag_asc');
while strcmp(happy_var,'no')
    centers_desc = [];
    imshow(im)
    if ~isempty(centers_asc)
        viscircles(centers_asc, radii_asc,'Color','r');
    end
    set(gcf,'windowbuttondownfcn',{@starttrack});
    waitfor(0,'userdata')
    desc_boundaries=get(0,'userdata');
    x_desc_boundaries = desc_boundaries(:,1); 
    y_desc_boundaries = handles.axes1.YLim(2)-desc_boundaries(:,2); 
    [xo yo radii_desc] = circle_fit(x_desc_boundaries,y_desc_boundaries);
    viscircles([xo yo], radii_desc,'Color','b'), viscircles(centers_asc, radii_asc,'Color','r');
    
    set(handles.Happy_str, 'Visible', 'on');set(handles.Yes_button, 'Visible', 'on');set(handles.No_button, 'Visible', 'on');
    waitforbuttonpress;
    pause(.4)
    if strcmp(happy_var,'yes')
        happy_desc = 'yes';
        centers_desc = [xo yo];
        set(handles.Happy_str, 'Visible', 'off');set(handles.Yes_button, 'Visible', 'off');set(handles.No_button, 'Visible', 'off');
        if strcmp(happy_asc,'yes')
            uiresume
        elseif flag_asc ==0
            centers_asc = [0 0]; radii_asc = 0;
            uiresume
        else
            break
        end
    else
        set(handles.Happy_str, 'Visible', 'off');set(handles.Yes_button, 'Visible', 'off');set(handles.No_button, 'Visible', 'off');
    end
end

%- Store variable in handle
guidata(hObject, handles);

% --- Executes on button press in Yes_button.
function Yes_button_Callback(hObject, eventdata, handles)
global happy_var
happy_var = 'yes';

% --- Executes on button press in No_button.
function No_button_Callback(hObject, eventdata, handles)
global happy_var
happy_var = 'no';

% --- Executes on button press in Forgo_Asc_Aorta_button.
function Forgo_Asc_Aorta_button_Callback(hObject, eventdata, handles)
global happy_asc happy_desc centers_asc radii_asc
centers_asc = [0 0]; radii_asc = 0;
happy_asc = 'yes';
setappdata(0,'happy_asc','yes')
if strcmp(happy_desc,'yes')
    uiresume
else
    return
end

% --- Executes on button press in Move_on_button.
function Move_on_button_Callback(hObject, eventdata, handles)
% hObject    handle to Move_on_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global flag_quit
flag_quit = 1;
uiresume
%- Store variable in handle
% guidata(hObject, handles);

% --- Executes on button press in Fuse_parasagittal_scan.
function Fuse_parasagittal_scan_Callback(hObject, eventdata, handles)
% hObject    handle to Fuse_parasagittal_scan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global flag_quit
flag_quit = 2;
GUI_PWV_parasagittal_fuse
uiresume
%- Store variable in handle
% guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- ADDITIONAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [centers_desc radii_desc centers_asc radii_asc] = find_desc_asc_from_ref_scans(im,centers_desc_ref,radii_desc_ref,centers_asc_ref,radii_asc_ref);
global ObjectPolarity_mode
    I = im;
    [centers, radii] = imfindcircles(I,[1 1000],'ObjectPolarity',ObjectPolarity_mode);
    %- Delete too large radii
    ind = find(radii>15);
    centers(ind,:) = [];
    radii(ind) = [];

    %% Descending aorta
    if isempty(radii_desc_ref)==1
        %- Y limits of ROI
        ymin_desc = round(size(I,1)*0.5); ymax_desc = round(size(I,1)*1);
        %- X limits of ROI
        xmin_desc = round(size(I,2)*0.4); xmax_desc = round(size(I,2)*0.6);
        %- Find circles within ROI
        centers_desc_potential = centers;radii_desc_potential = radii;
        for k=size(centers_desc_potential,1):-1:1
            if centers_desc_potential(k,1)<xmin_desc || centers_desc_potential(k,1)>xmax_desc || centers_desc_potential(k,2)<ymin_desc || centers_desc_potential(k,2)>ymax_desc
                centers_desc_potential(k,:) = [];
                radii_desc_potential(k,:) = [];
            end
        end
%     imshow(I),rectangle('Position',[xmin_desc ymin_desc xmax_desc-xmin_desc ymax_desc-ymin_desc]),viscircles(centers_desc_potential,radii_desc_potential,'Color','b')
        %- Select descending aorta
        [radii_desc imax] = max(radii_desc_potential);
        centers_desc = centers_desc_potential(imax,:);
    else
        diff_centers = []; diff_radii = []; diff_total = [];  
        for k=1:length(radii)
            diff_centers(k) = sqrt((centers(k,1)-centers_desc_ref(1))^2+(centers(k,2)-centers_desc_ref(2))^2);
            diff_radii(k) = radii(k) - radii_desc_ref;
            diff_total(k) = sqrt(diff_centers(k)^2 + diff_radii(k)^2);
        end
        [c imin] = min(diff_total);
        radii_desc = radii(imin);
        centers_desc = centers(imin,:)
    end
    
    %% Ascending aorta
    if isempty(radii_asc_ref)==1
        %- Y limits of ROI
        ymin_asc = 1; ymax_asc = round(size(I,1)*0.4);
        %- X limits of ROI
        xmin_asc = round(size(I,2)*0.25); xmax_asc = round(size(I,2)*0.5);
        %- Find circles within ROI
        centers_asc_potential = centers;radii_asc_potential = radii;
        for k=size(centers_asc_potential,1):-1:1
            if centers_asc_potential(k,1)<xmin_asc || centers_asc_potential(k,1)>xmax_asc || centers_asc_potential(k,2)<ymin_asc || centers_asc_potential(k,2)>ymax_asc
                centers_asc_potential(k,:) = [];
                radii_asc_potential(k,:) = [];
            end
        end
%     imshow(I),rectangle('Position',[xmin_asc ymin_asc xmax_asc-xmin_asc ymax_asc-ymin_asc]),viscircles(centers_asc_potential,radii_asc_potential,'Color','b')
        %- Select descending aorta
        [radii_asc imax] = max(radii_asc_potential);
        centers_asc = centers_asc_potential(imax,:);
    else
        diff_centers = []; diff_radii = []; diff_total = [];
        for k=1:length(radii)
            diff_centers(k) = sqrt((centers(k,1)-centers_asc_ref(1))^2+(centers(k,2)-centers_asc_ref(2))^2);
            diff_radii(k) = radii(k) - radii_asc_ref;
            diff_total(k) = diff_centers(k) + diff_radii(k);
        end
        [c imin] = min(diff_total);
        radii_asc = radii(imin);
        centers_asc = centers(imin,:)
    end
    
    
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

function starttrack(imagefig,varargins)
% global axsize
% imagefig.Units = 'pixels';
% imagefig.Position = [imagefig.Position(1) axsize(2) axsize(3) axsize(4)]
% get(0,'PointerLocation')
disp('tracking started')
set(gcf,'windowbuttondownfcn',{@stoptrack},'userdata',[]);
set(gcf,'windowbuttonmotionfcn',...
             'set(gcf,''userdata'',[get(gcf,''userdata'');get(0,''pointerlocation'')])');

function stoptrack(imagefig,varargins)

set(gcf,'windowbuttonmotionfcn',[]);

disp('tracking stopped')

units0=get(0,'units');
unitsf=get(gcf,'units');
unitsa=get(gca,'units');
                           
set(0,'units','pixels');
set(gcf,'units','pixels');
set(gca,'units','pixels');

x=get(gca,'xlim');
y=get(gca,'ylim');
axsize = plotboxpos(gca);
figsize=get(gcf,'position'); 
ratio=diag([diff(x)/axsize(3) diff(y)/axsize(4)]);
shift=figsize(1:2)+[axsize(1) axsize(2)];

set(0,'units',units0);
set(gcf,'units',unitsf);
set(gca,'units',unitsa);

data = get(gcf,'userdata');
mousetrail=(data-repmat(shift,size(data,1),1))*ratio;

set(gcf,'windowbuttondownfcn',{@starttrack});
set(0,'userdata',mousetrail);

function pos = plotboxpos(h)
%PLOTBOXPOS Returns the position of the plotted axis region
%
% pos = plotboxpos(h)
%
% This function returns the position of the plotted region of an axis,
% which may differ from the actual axis position, depending on the axis
% limits, data aspect ratio, and plot box aspect ratio.  The position is
% returned in the same units as the those used to define the axis itself.
% This function can only be used for a 2D plot.  
%
% Input variables:
%
%   h:      axis handle of a 2D axis (if ommitted, current axis is used).
%
% Output variables:
%
%   pos:    four-element position vector, in same units as h

% Copyright 2010 Kelly Kearney

% Check input

if nargin < 1
    h = gca;
end

if ~ishandle(h) || ~strcmp(get(h,'type'), 'axes')
    error('Input must be an axis handle');
end

% Get position of axis in pixels

currunit = get(h, 'units');
axisPos  = getpixelposition(h);

% Calculate box position based axis limits and aspect ratios

darismanual  = strcmpi(get(h, 'DataAspectRatioMode'),    'manual');
pbarismanual = strcmpi(get(h, 'PlotBoxAspectRatioMode'), 'manual');

if ~darismanual && ~pbarismanual
    
    pos = axisPos;
    
else

    xlim = get(h, 'XLim');
    ylim = get(h, 'YLim');
    
    % Deal with axis asclimits auto-set via Inf/-Inf use
    
    if any(isinf([xlim ylim]))
        hc = get(h, 'Children');
        hc(~arrayfun( @(h) isprop(h, 'XData' ) & isprop(h, 'YData' ), hc)) = [];
        xdata = get(hc, 'XData');
        if iscell(xdata)
            xdata = cellfun(@(x) x(:), xdata, 'uni', 0);
            xdata = cat(1, xdata{:});
        end
        ydata = get(hc, 'YData');
        if iscell(ydata)
            ydata = cellfun(@(x) x(:), ydata, 'uni', 0);
            ydata = cat(1, ydata{:});
        end
        isplotted = ~isinf(xdata) & ~isnan(xdata) & ...
                    ~isinf(ydata) & ~isnan(ydata);
        xdata = xdata(isplotted);
        ydata = ydata(isplotted);
        if isempty(xdata)
            xdata = [0 1];
        end
        if isempty(ydata)
            ydata = [0 1];
        end
        if isinf(xlim(1))
            xlim(1) = min(xdata);
        end
        if isinf(xlim(2))
            xlim(2) = max(xdata);
        end
        if isinf(ylim(1))
            ylim(1) = min(ydata);
        end
        if isinf(ylim(2))
            ylim(2) = max(ydata);
        end
    end

    dx = diff(xlim);
    dy = diff(ylim);
    dar = get(h, 'DataAspectRatio');
    pbar = get(h, 'PlotBoxAspectRatio');

    limDarRatio = (dx/dar(1))/(dy/dar(2));
    pbarRatio = pbar(1)/pbar(2);
    axisRatio = axisPos(3)/axisPos(4);

    if darismanual
        if limDarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/limDarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * limDarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    elseif pbarismanual
        if pbarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/pbarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * pbarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    end
end

% Convert plot box position to the units used by the axis

hparent = get(h, 'parent');
hfig = ancestor(hparent, 'figure'); % in case in panel or similar
currax = get(hfig, 'currentaxes');

temp = axes('Units', 'Pixels', 'Position', pos, 'Visible', 'off', 'parent', hparent);
set(temp, 'Units', currunit);
pos = get(temp, 'position');
delete(temp);

set(hfig, 'currentaxes', currax);

function [xo yo R] = circle_fit(x,y)
% A function to find the best circle fit (radius and center location) to
% given x,y pairs
% 
% Val Schmidt
% Center for Coastal and Ocean Mapping
% University of New Hampshire
% 2012
% 
% Arguments:
% x:         x coordinates
% y:         y coordinates
%
% Output:
% xo:        circle x coordinate center
% yo:        circle y coordinate center
% R:         circle radius
x = x(:);
y = y(:);
% Fitting a circle to the data - least squares style. 
%Here we start with
% (x-xo).^2 + (y-yo).^2 = R.^2
% Rewrite:
% x.^2 -2 x xo + xo^2 + y.^2 -2 y yo + yo.^2 = R.^2
% Put in matrix form:
% [-2x -2y 1 ] [xo yo -R^2+xo^2+yo^2]' = -(x.^2 + y.^2)
% Solve in least squares way...
A = [-2*x -2*y ones(length(x),1)];
x = A\-(x.^2+y.^2);
xo=x(1);
yo=x(2);
R = sqrt(  xo.^2 + yo.^2  - x(3));

