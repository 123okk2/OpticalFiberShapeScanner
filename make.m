function varargout = make(varargin)
% MAKE MATLAB code for make.fig
%      MAKE, by itself, creates a new MAKE or raises the existing
%      singleton*.
%
%      H = MAKE returns the handle to a new MAKE or the handle to
%      the existing singleton*.
%
%      MAKE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAKE.M with the given input arguments.
%
%      MAKE('Property','Value',...) creates a new MAKE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before make_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to make_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help make

% Last Modified by GUIDE v2.5 04-Feb-2015 10:24:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @make_OpeningFcn, ...
                   'gui_OutputFcn',  @make_OutputFcn, ...
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


% --- Executes just before make is made visible.
function make_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to make (see VARARGIN)

% Choose default command line output for make
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes make wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = make_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in play.
function play_Callback(hObject, eventdata, handles)
% hObject    handle to play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
colordef black
global iter
global status
global WLref
global WLin_a
global save

DN=str2double(get(handles.dn,'string'));
nodeFBG=str2double(get(handles.nodefbg,'string'));
h1=str2double(get(handles.mid,'string'));
h=ones(1,nodeFBG+1).*h1;
h(1) = str2double(get(handles.left,'string'));
h(nodeFBG+1) = str2double(get(handles.right,'string'));
FileName = get(handles.filename, 'string');
WLref = csvread('ref_0001.csv', 531, 0, [531, 0, 530+nodeFBG*3, 0]);
d =[0.0796	0.0783	0.0804;
    0.0769	0.0863	0.0732;
    0.0749	0.0796	0.0605;
    0.0959	0.0619	0.0761];

ang = [	2.5333    0.5261   -1.6939;
        2.7369    0.7555   -1.5296;
        3.0208    1.0560   -1.2374;
        -2.7948    1.5246   -1.0304];

d_A=d(:,1); d_B=d(:,2); d_C=d(:,3);
ang_A=ang(:,1); ang_B=ang(:,2); ang_C=ang(:,3);
if status~=0
    iter=1;
else
    iter=save;
end
status=1;

s = [linspace(0, h(1), DN)];
disp(s);
for i=2:nodeFBG+1
    s = [s linspace(sum(h(1:i-1))+h(i)/DN, sum(h(1:i)), DN)];
end

while iter ~=0
count = 0;
while (exist(sprintf('%s_%04d.csv', FileName, iter),'file') ~= 2)
    count = count + 1;
    pause(0.01);
    if count == 500
        iter= -1;
        break
    end
end
if iter ~=-1
   PeakN= csvread(sprintf('%s_%04d.csv', FileName, iter), 529, 1, [529, 1, 529, 1]) ;
   WLin = csvread(sprintf('%s_%04d.csv', FileName, iter), 531, 0, [531, 0, 530+PeakN, 0]);
   while PeakN>12
        [DPWL DP]=min(WLin(2:PeakN)-WLin(1:PeakN-1));
        WLin=[WLin(1:DP-1)' mean(WLin(DP:DP+1)) WLin(DP+2:PeakN)']';
        PeakN=PeakN-1;
    end
end

if PeakN==12
strain=(WLin./WLref - 1);
k_A = strain(1:4)./d_A;
k_B = strain(5:8)./d_B;
k_C = strain(9:12)./d_C;
for i=1:4
    eq1(i,:)=[sin(ang_A(i)) cos(ang_A(i)); sin(ang_B(i)) cos(ang_B(i))]^-1*[k_A(i); k_B(i)];
    eq2(i,:)=[sin(ang_A(i)) cos(ang_A(i)); sin(ang_C(i)) cos(ang_C(i))]^-1*[k_A(i); k_C(i)];
    eq3(i,:)=[sin(ang_B(i)) cos(ang_B(i)); sin(ang_C(i)) cos(ang_C(i))]^-1*[k_B(i); k_C(i)];
    
    [TH1(i),K1(i)]=cart2pol(eq1(i,1),eq1(i,2));
    [TH2(i),K2(i)]=cart2pol(eq2(i,1),eq2(i,2));
    [TH3(i),K3(i)]=cart2pol(eq3(i,1),eq3(i,2));
end
eq=(eq1+eq2+eq3)/3;
[THin,kin]=cart2pol(eq(:,1),eq(:,2));
%%%%%%%%%%%%%%%%    ShapeReconstruction_Circle %%%%%%%%%%%%%%%%%%%%
R=ShapeReconst_circle(DN,nodeFBG,kin,THin,h);
r(:,1)=mean(R(:,:,1)')';
r(:,2)=mean(R(:,:,2)')';
r(:,3)=mean(R(:,:,3)')';
tip_direction(1,:)=(r((nodeFBG+1)*DN,:)-r((nodeFBG+1)*DN-1,:))/norm((r((nodeFBG+1)*DN,:)-r((nodeFBG+1)*DN-1,:)));
%%%%%%%%%%%%%%%%%%%%%% Surf plot %%%%%%%%%%%%%%%%%%%%%%%%%%

if status==1

    Rgraph=surf(handles.axes1,R(:,:,3),R(:,:,1),R(:,:,2),'edgecolor','none');
    set(Rgraph,'LineStyle','none','FaceColor',[0 1 0]);
    camlight (156,110);
    lighting GOURAUD
    material shiny
    
    x=get(handles.slider3,'value');
    y=get(handles.slider4,'value');
    
    hold on
    if y>=0
        if (x>90 && x<270)
            plot3(handles.axes1,zeros(1,500),r(:,1),r(:,2),'w','linewidth',1);
            plot3(handles.axes1,r(:,3),zeros(1,500)-50,r(:,2),'w','linewidth',1);
            plot3(handles.axes1,r(:,3),r(:,1),zeros(1,500)-50,'w','linewidth',1);
        else
            plot3(handles.axes1,zeros(1,500),r(:,1),r(:,2),'w','linewidth',1);
            plot3(handles.axes1,r(:,3),zeros(1,500)+50,r(:,2),'w','linewidth',1);
            plot3(handles.axes1,r(:,3),r(:,1),zeros(1,500)-50,'w','linewidth',1);
    end
    else
        if (x<90 || x>270)
            plot3(handles.axes1,zeros(1,500),r(:,1),r(:,2),'w','linewidth',1);
            plot3(handles.axes1,r(:,3),zeros(1,500)+50,r(:,2),'w','linewidth',1);
            plot3(handles.axes1,r(:,3),r(:,1),zeros(1,500)+50,'w','linewidth',1);
        else
            plot3(handles.axes1,zeros(1,500),r(:,1),r(:,2),'w','linewidth',1);
            plot3(handles.axes1,r(:,3),zeros(1,500)-50,r(:,2),'w','linewidth',1);
            plot3(handles.axes1,r(:,3),r(:,1),zeros(1,500)+50,'w','linewidth',1);
        end
    end
    set(handles.tippoint,'String',sprintf('\t%0.2f\t%0.2f\t%0.2f',r((nodeFBG+1)*DN,1),r((nodeFBG+1)*DN,2),r((nodeFBG+1)*DN,3)));
    set(handles.tipdirection,'String',sprintf('\t%0.2f\t%0.2f\t%0.2f',tip_direction(1,1),tip_direction(1,2),tip_direction(1,3)));
    
    hold off
    set(handles.az,'string',sprintf('%0.2f',x));
    set(handles.el,'string',sprintf('%0.1f',y));
    set(handles.axes1,'xlim',[0 115],'ylim',[-50 50],'zlim',[-50 50],'view',[x y]);
    drawnow
end
end
iter = iter+1
end


% --- Executes on button press in changeref.
function changeref_Callback(hObject, eventdata, handles)
% hObject    handle to changeref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global WLref
global WLin_A
WLref=WLin_A;

% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global iter
global status
iter=-1;
status=-1;

% --- Executes on button press in pause.
function pause_Callback(hObject, eventdata, handles)
% hObject    handle to pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global iter
global status
global save
save=iter;
iter=-1;
status=0;


function filename_Callback(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filename as text
%        str2double(get(hObject,'String')) returns contents of filename as a double


% --- Executes during object creation, after setting all properties.
function filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nodefbg_Callback(hObject, eventdata, handles)
% hObject    handle to nodefbg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nodefbg as text
%        str2double(get(hObject,'String')) returns contents of nodefbg as a double


% --- Executes during object creation, after setting all properties.
function nodefbg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nodefbg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dn_Callback(hObject, eventdata, handles)
% hObject    handle to dn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dn as text
%        str2double(get(hObject,'String')) returns contents of dn as a double


% --- Executes during object creation, after setting all properties.
function dn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function left_Callback(hObject, eventdata, handles)
% hObject    handle to left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of left as text
%        str2double(get(hObject,'String')) returns contents of left as a double


% --- Executes during object creation, after setting all properties.
function left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mid_Callback(hObject, eventdata, handles)
% hObject    handle to mid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mid as text
%        str2double(get(hObject,'String')) returns contents of mid as a double


% --- Executes during object creation, after setting all properties.
function mid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function right_Callback(hObject, eventdata, handles)
% hObject    handle to right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of right as text
%        str2double(get(hObject,'String')) returns contents of right as a double


% --- Executes during object creation, after setting all properties.
function right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tippoint_Callback(hObject, eventdata, handles)
% hObject    handle to tippoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tippoint as text
%        str2double(get(hObject,'String')) returns contents of tippoint as a double


% --- Executes during object creation, after setting all properties.
function tippoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tippoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tipdirection_Callback(hObject, eventdata, handles)
% hObject    handle to tipdirection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tipdirection as text
%        str2double(get(hObject,'String')) returns contents of tipdirection as a double

% --- Executes during object creation, after setting all properties.
function tipdirection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tipdirection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x_Callback(hObject, eventdata, handles)
% hObject    handle to x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x as text
%        str2double(get(hObject,'String')) returns contents of x as a double


% --- Executes during object creation, after setting all properties.
function x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function y_Callback(hObject, eventdata, handles)
% hObject    handle to y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y as text
%        str2double(get(hObject,'String')) returns contents of y as a double


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function az_Callback(hObject, eventdata, handles)
% hObject    handle to az (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of az as text
%        str2double(get(hObject,'String')) returns contents of az as a double


% --- Executes during object creation, after setting all properties.
function az_CreateFcn(hObject, eventdata, handles)
% hObject    handle to az (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function el_Callback(hObject, eventdata, handles)
% hObject    handle to el (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of el as text
%        str2double(get(hObject,'String')) returns contents of el as a double


% --- Executes during object creation, after setting all properties.
function el_CreateFcn(hObject, eventdata, handles)
% hObject    handle to el (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
