function varargout = GUI_ref_top_aorta(varargin)
% GUI_REF_TOP_AORTA MATLAB code for GUI_ref_top_aorta.fig
%      GUI_REF_TOP_AORTA, by itself, creates a new GUI_REF_TOP_AORTA or raises the existing
%      singleton*.
%
%      H = GUI_REF_TOP_AORTA returns the handle to a new GUI_REF_TOP_AORTA or the handle to
%      the existing singleton*.
%
%      GUI_REF_TOP_AORTA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_REF_TOP_AORTA.M with the given input arguments.
%
%      GUI_REF_TOP_AORTA('Property','Value',...) creates a new GUI_REF_TOP_AORTA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_ref_top_aorta_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_ref_top_aorta_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_ref_top_aorta

% Last Modified by GUIDE v2.5 08-Dec-2020 14:46:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_ref_top_aorta_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_ref_top_aorta_OutputFcn, ...
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


% --- Executes just before GUI_ref_top_aorta is made visible.
function GUI_ref_top_aorta_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_ref_top_aorta (see VARARGIN)
global happy_var happy_aorta Im_top centers_aorta happy_automatic ObjectPolarity_mode flag_quit
happy_var = 'no'; happy_aorta = 'no'; happy_automatic = 'no'; centers_aorta = []; ObjectPolarity_mode = 'bright'; flag_quit = 0;
handles.output = hObject;
Im_top = getappdata(0,'Im_top');
axes(handles.axes1), imshow(Im_top)
set(handles.Happy_str, 'Visible', 'off');set(handles.Button_yes, 'Visible', 'off');set(handles.Button_no, 'Visible', 'off');

% Choose default command line output for GUI_ref_asc_desc
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_ref_asc_desc wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = GUI_ref_top_aorta_OutputFcn(hObject, eventdata, handles) 
global centers radii happy_var flag_quit
setappdata(0,'centers_top',centers)
setappdata(0,'radii_top',radii)

if flag_quit==1
    setappdata(0,'flag_quit',1)
    uiresume
else
    setappdata(0,'flag_quit',0)
end

if strcmp(happy_var,'yes')
    uiresume
    close
end


varargout{1} = handles.output;

% --- Executes on button press in Button_automatic.
function Button_automatic_Callback(hObject, eventdata, handles)
global Im_top ObjectPolarity_mode happy_var centers radii
I = Im_top(round(size(Im_top,1)*.6):round(size(Im_top,1)*1),round(size(Im_top,2)*.3):round(size(Im_top,2)*.7));
[centers, radii] = imfindcircles(I,[3 30],'ObjectPolarity',ObjectPolarity_mode);
%- Extract ref aorta and radius
[c imax] = max(radii);
centers = centers(imax,:);
centers = [centers(1)+round(size(Im_top,2)*.3) centers(2)+round(size(Im_top,1)*.6)]
radii = radii(imax);
happy_var = 'no';
if strcmp(happy_var,'no')
    axes(handles.axes1),imshow(Im_top), viscircles(centers, radii,'Color','r')
    set(handles.Happy_str, 'Visible', 'on');set(handles.Button_yes, 'Visible', 'on');set(handles.Button_no, 'Visible', 'on');
    waitforbuttonpress;
    pause(.4)
    if strcmp(happy_var,'yes')
        uiresume
    else
        pause(.4)
        centers = [], radii = [];
        imshow(Im_top)
        set(handles.Happy_str, 'Visible', 'off');set(handles.Button_yes, 'Visible', 'off');set(handles.Button_no, 'Visible', 'off');
    end
end

%- Store variable in handle
guidata(hObject, handles);


% --- Executes on button press in Button_manual.
function Button_manual_Callback(hObject, eventdata, handles)
global Im_top ObjectPolarity_mode happy_var centers radii
[centers, radii] = imfindcircles(Im_top,[5 30],'ObjectPolarity',ObjectPolarity_mode);
%- Delete too large radii
ind1 = find(radii>20); ind2 = find(radii<3); ind = [ind1 ; ind2]
centers(ind,:) = [];
radii(ind) = [];
%% Ascending aorta
happy_var = 'no';
while strcmp(happy_var,'no')
    imshow(Im_top)
    set(gcf,'windowbuttondownfcn',{@starttrack});
    waitfor(0,'userdata')
    aorta_boundaries=get(0,'userdata');
    x_aorta_boundaries = aorta_boundaries(:,1); 
    y_aorta_boundaries = handles.axes1.YLim(2)-aorta_boundaries(:,2); 
    [xo yo radii_aorta] = circle_fit(x_aorta_boundaries,y_aorta_boundaries);
    viscircles([xo yo], radii_aorta,'Color','r')
    
    
%     imshow(Im_top), viscircles(centers, radii,'Color','r');
%     [x,y] = ginput(1);
%     [x_ref y_ref radii_ref] = find_ref(centers,radii,x,y);
%     imshow(Im_top), viscircles([x_ref y_ref], radii_ref,'Color','r')
    set(handles.Happy_str, 'Visible', 'on');set(handles.Button_yes, 'Visible', 'on');set(handles.Button_no, 'Visible', 'on');
    waitforbuttonpress;
    pause(.4)
    if strcmp(happy_var,'yes')
        centers = [xo yo]; radii = radii_aorta;
        happy_var = 'yes';
        set(handles.Happy_str, 'Visible', 'off');set(handles.Button_yes, 'Visible', 'off');set(handles.Button_no, 'Visible', 'off');
        uiresume
    else
        set(handles.Happy_str, 'Visible', 'off');set(handles.Button_yes, 'Visible', 'off');set(handles.Button_no, 'Visible', 'off');
    end
end

%- Store variable in handle
guidata(hObject, handles);

% --- Executes on button press in Button_yes.
function Button_yes_Callback(hObject, eventdata, handles)
global happy_var
happy_var = 'yes';

% --- Executes on button press in Button_no.
function Button_no_Callback(hObject, eventdata, handles)
global happy_var
happy_var = 'no';

% --- Executes on button press in Move_On_Button.
function Move_On_Button_Callback(hObject, eventdata, handles)
global flag_quit
flag_quit = 1;
uiresume
%- Store variable in handle
guidata(hObject, handles);


%%
%- ADDITIONAL FUNCTIONS
%%
function [x_ref y_ref rad_ref] = find_ref(centers,radii,x,y);
for k=1:size(centers)
   d(k) = sqrt((centers(k,1)-x)^2 + (centers(k,2)-y)^2);
end
[c imin] = min(d)
x_ref = centers(imin,1);
y_ref = centers(imin,2);
rad_ref = radii(imin);

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


