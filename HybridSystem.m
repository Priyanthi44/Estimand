function varargout = HybridSystem(varargin)
% HYBRIDSYSTEM MATLAB code for HybridSystem.fig
%      HYBRIDSYSTEM, by itself, creates a new HYBRIDSYSTEM or raises the existing
%      singleton*.
%
%      H = HYBRIDSYSTEM returns the handle to a new HYBRIDSYSTEM or the handle to
%      the existing singleton*.
%
%      HYBRIDSYSTEM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HYBRIDSYSTEM.M with the given input arguments.
%
%      HYBRIDSYSTEM('Property','Value',...) creates a new HYBRIDSYSTEM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HybridSystem_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HybridSystem_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HybridSystem

% Last Modified by GUIDE v2.5 12-Apr-2018 17:16:26

% Begin initialization code - DO NOT EDIT
size(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HybridSystem_OpeningFcn, ...
                   'gui_OutputFcn',  @HybridSystem_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
    %create a figure to house the GUI
% figure
% ellipse_position = [0.3 0.6 0.1 0.2];
% ellipse_h = annotation('ellipse',ellipse_position,...
%                     'facecolor', [1 0 0]);
%                  
%create an annotation object 
% ellipse_position = [0.4 0.6 0.1 0.2];
% ellipse_h = annotation('ellipse',ellipse_position,...
%                     'facecolor', [1 0 0]);
%                  
% %create an editable textbox object
% edit_box_h = uicontrol('style','edit',...
%                     'units', 'normalized',...
%                     'position', [0.3 0.4 0.4 0.1]);
%  
% %create a "push button" user interface (UI) object
% but_h = uicontrol('style', 'pushbutton',...
%                     'string', 'Update Color',...
%                     'units', 'normalized',...
%                     'position', [0.3 0 0.4 0.2],...
%                     'callback', {@eg_fun,edit_box_h, ellipse_h });
%              
% %Slider object to control ellipse size
% uicontrol('style','Slider',...        
%             'Min',0.5,'Max',2,'Value',1,...
%             'units','normalized',...
%             'position',[0.1    0.2    0.08    0.25],...
%             'callback',{@change_size,ellipse_h,ellipse_position });
%          
% uicontrol('Style','text',...
%             'units','normalized',...
%             'position',[0    0.45    0.2    0.1],...
%             'String','Ellipse Size')
end
% End initialization code - DO NOT EDIT


% --- Executes just before HybridSystem is made visible.
function HybridSystem_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HybridSystem (see VARARGIN)

% Choose default command line output for HybridSystem
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes HybridSystem wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HybridSystem_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% spikingneuralnetwork();
   
ellipse_position = [0.3 0.6 0.1 0.2];
ellipse_h = annotation('ellipse',ellipse_position,...
                    'facecolor', [1 0 0]);
                          

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(ellipse_h);
ellipse_position = [0.5 0.6 0.1 0.2];
ellipse_h = annotation('ellipse',ellipse_position,...
                    'facecolor', [1 0 0]);
                delete(ellipse_h);
                

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ellipse_position = [0.3 0.6 0.1 0.2];
ellipse_h = annotation('ellipse',ellipse_position,...
                    'facecolor', [1 0 0]);
