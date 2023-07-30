%TUMORSIM: A MATLAB program for simulation of tumor growth before and after
%irradiation
%created by Yoichi Watanabe, Ph.D. in October 2013.
%Version 5.4 last updated on 08/07/2015
%12/03/2013: The first version was completed.
%12/05/2013: add Cell clearance time
%-- assume vascular mass is lamda times tumor mass
%-- added parameer optimization function
%-- fun_vdp is in this m file.
%12/06/2013: add instant cell kill model option (use survival fraction)
%12/10/2013: 
%           Radiation-response model Model A: instant cell kill
%           Radiation-response model Model B: continous cell kill
%           Model A: the tumor volume = only proliferating cells
%           Model B: the tumor volume = proliferating + non-dividing cells
%           plot the clinical data on both axes 1 and axes 2
%01/16/2014: revise/correct minor errors
%04/23/2014: the first term in eq.(1) for y1 was revised by removing c(1)
%05/09/2014: change lambda equation to allow the volume change after
%           irradiation. Function vdp9a
%05/26/2014: 
%   assume that the rammbda_t is constant while the radiation is on.
%   added a function total_cellkill 
%10/31/2014: redefine the tumor doubling time: new_T = old_T /2
%11/02/2014: print the lamda on the day of GKSRS
%            odeset tolerance values were changed.
%08/06/2015: bugs with Model A were fixed.
%12/13/2015: replace one theta to (theta1 and theta2)
function varargout = TUMORSIM55(varargin)
% TUMORSIM55 MATLAB code for TUMORSIM55.fig
%      TUMORSIM55, by itself, creates a new TUMORSIM55 or raises the existing
%      singleton*.
%
%      H = TUMORSIM55 returns the handle to a new TUMORSIM55 or the handle to
%      the existing singleton*.
%
%      TUMORSIM55('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TUMORSIM55.M with the given input arguments.
%
%      TUMORSIM55('Property','Value',...) creates a new TUMORSIM55 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TUMORSIM55_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TUMORSIM55_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TUMORSIM55

% Last Modified by GUIDE v2.5 13-Dec-2015 15:35:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TUMORSIM55_OpeningFcn, ...
                   'gui_OutputFcn',  @TUMORSIM55_OutputFcn, ...
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


% --- Executes just before TUMORSIM55 is made visible.
function TUMORSIM55_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TUMORSIM55 (see VARARGIN)

% Choose default command line output for TUMORSIM55
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TUMORSIM55 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TUMORSIM55_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function alpha_Callback(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha as text
%        str2double(get(hObject,'String')) returns contents of alpha as a double


% --- Executes during object creation, after setting all properties.
function alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alphabeta_Callback(hObject, eventdata, handles)
% hObject    handle to alphabeta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alphabeta as text
%        str2double(get(hObject,'String')) returns contents of alphabeta as a double


% --- Executes during object creation, after setting all properties.
function alphabeta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphabeta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dose_Callback(hObject, eventdata, handles)
% hObject    handle to dose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dose as text
%        str2double(get(hObject,'String')) returns contents of dose as a double


% --- Executes during object creation, after setting all properties.
function dose_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rad_on_Callback(hObject, eventdata, handles)
% hObject    handle to rad_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rad_on as text
%        str2double(get(hObject,'String')) returns contents of rad_on as a double


% --- Executes during object creation, after setting all properties.
function rad_on_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rad_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_double_Callback(hObject, eventdata, handles)
% hObject    handle to T_double (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_double as text
%        str2double(get(hObject,'String')) returns contents of T_double as a double


% --- Executes during object creation, after setting all properties.
function T_double_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_double (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vasc_ret_factor_Callback(hObject, eventdata, handles)
% hObject    handle to vasc_ret_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vasc_ret_factor as text
%        str2double(get(hObject,'String')) returns contents of vasc_ret_factor as a double


% --- Executes during object creation, after setting all properties.
function vasc_ret_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vasc_ret_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vasc_tumor_ratio_Callback(hObject, eventdata, handles)
% hObject    handle to vasc_tumor_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vasc_tumor_ratio as text
%        str2double(get(hObject,'String')) returns contents of vasc_tumor_ratio as a double


% --- Executes during object creation, after setting all properties.
function vasc_tumor_ratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vasc_tumor_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_step_Callback(hObject, eventdata, handles)
% hObject    handle to T_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_step as text
%        str2double(get(hObject,'String')) returns contents of T_step as a double


% --- Executes during object creation, after setting all properties.
function T_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_duration_Callback(hObject, eventdata, handles)
% hObject    handle to T_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_duration as text
%        str2double(get(hObject,'String')) returns contents of T_duration as a double


% --- Executes during object creation, after setting all properties.
function T_duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_cycle_Callback(hObject, eventdata, handles)
% hObject    handle to T_cycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_cycle as text
%        str2double(get(hObject,'String')) returns contents of T_cycle as a double


% --- Executes during object creation, after setting all properties.
function T_cycle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_cycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function n_cycle_Callback(hObject, eventdata, handles)
% hObject    handle to n_cycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n_cycle as text
%        str2double(get(hObject,'String')) returns contents of n_cycle as a double


% --- Executes during object creation, after setting all properties.
function n_cycle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_cycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tag_start.
function tag_start_Callback(hObject, eventdata, handles)
% hObject    handle to tag_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%point model of tumor growth under irradiation
%yy1 = yy2 + yy3
%yy2 = mass of dividing cells
%yy3 = mass of non-dividing cells
%yy4 = lamda of tumor
%yy5 = mass of blood inside the tumor
%yy6 = effective lamda modulated by radiation effect
%yy7 = yy5/yy2

%d = dose [Gy]
%Time = real time/cell cycle time
%the cell cycle time T* (T_star)- 1 day
%global a b eta kappa mu theta dose p q rad_on expdatain voldata
global p q rad_on expdatain voldata eta_clearance instant_kill_model

%=================== input parameters ======
dose = str2num(get(handles.dose,'String')); %single dose [Gy]
T_double = str2num(get(handles.T_double,'String')); %days singular at alog2
T_repair = str2num(get(handles.T_repair,'String')); %repair time in days
T_star = str2num(get(handles.T_cycle,'String')); %day cell cycle time
T_clearance = str2num(get(handles.T_clear,'String')); %day cell clearance time
theta = str2num(get(handles.vasc_ret_factor,'String')); %relative ratio of lamda_v and lamda
alpha = str2num(get(handles.alpha,'String')); %1/Gy
alpha_over_beta = str2num(get(handles.alphabeta,'String'));
kappa = str2num(get(handles.kappa,'String'));
nn = str2num(get(handles.n_cycle,'String')); % the cell cyle number at which the colonies are counted
dt = str2num(get(handles.T_step,'String')); %time step for simulation [day]
T_duration = str2num(get(handles.T_duration,'String')); %the total time in days
tmax = T_duration(2);
Y0 = str2num(get(handles.Y_initial,'String'));
rad_effective = str2num(get(handles.rad_on,'String')); %hours when the radiation is killing the cells
rad_on = rad_effective*(1/T_star); %turn on radiation, single instant dose
pause on;
%===========================================
plot_exp_data_on = get(handles.pltexpdata,'Max');
plot_exp_data_option = get(handles.pltexpdata,'Value');
on_semilogY = str2num(get(handles.on_semilogy,'String'));
Instant_cell_kill_on = get(handles.instantkill,'Max');
Instant_cell_kill_option = get(handles.instantkill,'Value');
cont_cell_kill_on = get(handles.cont_cell_kill,'Max');
cont_cell_kill_option = get(handles.cont_cell_kill,'Value');

%ode solver tolerance updated on 11/2/2014
ode_options = odeset('RelTol',1e-9,'AbsTol',[1e-9 1e-9 1e-9]);

if Instant_cell_kill_option == Instant_cell_kill_on
   instant_kill_model = 1;
else
   instant_kill_model = 0; 
end
if cont_cell_kill_option == cont_cell_kill_on
   cont_kill_model = 1;
else
   cont_kill_model = 0; 
end

if cont_kill_model + instant_kill_model ~= 1
    msgbox('Choose Model A or B');
    pause;
    exit;
end

if expdatain ~= 1
expdatain = 0; %no experimentl data input (default)
end
beta = alpha/alpha_over_beta;
chi = alpha*dose+beta*dose*dose;
nn_min = chi/3.;
if nn >= nn_min
qsi = chi/(3*nn);
q = qsi;
else
    msgbox('nn is too small');
    nn_min
    pause;
    q = 0.0;
end
p = 1 - q;
alog2 = log(2);
ii = size(T_double,2);
for i=1:ii
lammda(i) = alog2/(2*T_double(i))*T_star; %lammda*T*
const_lammda(i) = lammda(i)/kappa;
end
lammda_v = lammda*theta(1); %set the lamgda at t=0
eta_clearance = alog2/T_clearance*T_star; %cell clearance rate

i1 = floor(tmax/dt);
for j=1:1:i1
  Tspan(j) = (dt/T_star)*(j-1);
end

%---------------Select model B or A --------------
if (cont_kill_model == 1)
%Model B
    for j=1:1:i1
        t = Tspan(j);
        if (t >= rad_on(1) && t <= rad_on(2))
            pp(j) = p;
        else
            pp(j) = 1.0;
        end
    end
    j1 = floor(rad_on(1)/dt);
    j2 = i1 - j1 + 1;
    for j = 1:1:j1
        TTspan1(j) = (dt/T_star)*(j-1);
    end
    for j = j1:1:i1
        TTspan2(j-j1+1) = (dt/T_star)*(j-1);
    end
elseif (instant_kill_model == 1)
%Model A
    for j=1:1:i1
        j1 = floor(rad_on(1)/dt);
        j2 = i1 - j1 + 1;
        for j = 1:1:j1
            TTspan1(j) = (dt/T_star)*(j-1);
            pp(j) = 1.0;
        end
        for j = j1:1:i1
            TTspan2(j-j1+1) = (dt/T_star)*(j-1);
            pp(j) = 1.0;
        end 
    end
end

if dose == 0 
    mu = 0; %if no irradiation, no repair
else
    mu = alog2/T_repair*T_star;
end

%---------------------------------
% y = mass/mass_nominal
% T = time in days/T_star
%Y(1) = yy2
%Y(2) = yy3
%Y(3) = yy4
%tt = zeros(i1,1);
%y = zeros(3,1);

yy = zeros(i1,7,ii);
%loop for ii of Td ---------------------------------------
for i=1:1:ii
c0(1) = lammda(i); %at the first diagnostic exam
c0(2) = mu;
c0(3) = theta(1);
Y0(3) = c0(1); %the initial condition of lambda(t)

if (cont_kill_model == 1)
    [T, Y] = ode15s(@(t,y)fun_vdp9a(t,y,c0),TTspan1,Y0, ode_options);
    YY0(1) = Y(j1,1);
    YY0(2) = 0.0;
    YY0(3) = Y(j1,3);
    c0(1) = Y(j1,3);
    c0(2) = mu;
    c0(3) = theta(2);
    [T1, Y1] = ode15s(@(t,y)fun_vdp9a(t,y,c0),TTspan2,YY0, ode_options);
    for j=j1+1:1:i1
        T(j) = T1(j-j1);
        Y(j,:) = Y1(j-j1,:);
    end
elseif (instant_kill_model == 1)
    [T, Y] = ode15s(@(t,y)fun_vdp9a(t,y,c0),TTspan1,Y0, ode_options);
    ss = exp(-chi);
    YY0(1) = ss*Y(j1,1);
    YY0(2) = (1-ss)*Y(j1,1);
    YY0(3) = Y(j1,3);
    [T1, Y1] = ode15s(@(t,y)fun_vdp9a(t,y,c0),TTspan2,YY0, ode_options);
    for j=j1+1:1:i1
        T(j) = T1(j-j1);
        Y(j,:) = Y1(j-j1,:);
    end
else
end

for k=1:1:2
    yy(:,k+1,i) = Y(:,k);
end
yy(:,4,i) = Y(:,3)/T_star;
for l=1:1:i1
    yy(l,5,i) = yy(l,4,i)*yy(l,2,i)/const_lammda(i);
    yy(l,6,i) = alog2/(pp(l)*yy(l,4,i));
    yy(l,7,i) = yy(l,5,i)/yy(l,2,i);
end
    yy(:,1,i) = yy(:,2,i) + yy(:,3,i); %total tumore volume

if (cont_kill_model == 1)
    total_cellkill2(i) = intgl_killterm(rad_on, T, yy(:,2,i), q, dt);
    total_cellkill(i) = yy(rad_on(1),2,i) - yy(rad_on(2),2,i);
    cellsurvival_rate(i) = 1-total_cellkill(i)/yy(rad_on(1),2,i);
end

end 
%---------- loop for ii (Td)-----------------------------------

% plot the results
plotparam = struct('ytitle',...
    {'Total tumor volume',...
    'Proliferating Tumor volume',...
    'Volume of non-dividing cells',...   
    'lamda',...
    'Vasculatur volume',...
    'Effective tumor doubling time [days]',...
    'Vasculatur volume/tumor volume'...
    });
lineopt = struct('linetype',{'-r','-b','-g'});

tm = T*T_star;
for no_plt = 1:1:2
    switch no_plt
        case 1
            axes(handles.axes1);
        case 2
            axes(handles.axes2);
    end
    cla   
    hold all;
    for i=1:1:ii
    plot(tm, yy(:,no_plt,i),lineopt(i).linetype,'LineWidth',2);
        if on_semilogY(no_plt) == 1
        set(gca,'yscale','log') %use for semilogy
        else
          set(gca,'yscale','linear')  
        end
    end
    if plot_exp_data_option == plot_exp_data_on && expdatain == 1
    plot(voldata(:,1), voldata(:,2),'dr','MarkerSize', 10, 'MarkerFaceColor','r');
        if on_semilogY(no_plt) == 1
        set(gca,'yscale','log') %use for semilogy
        else
          set(gca,'yscale','linear')  
        end
    end
    
    xlabel('Time [days]');ylabel(plotparam(no_plt).ytitle);
    if plot_exp_data_option == plot_exp_data_on && expdatain == 1
        hleg1=legend(num2str(T_double(1)),num2str(T_double(2)),...
            num2str(T_double(3)),'Exp.Data','Location','NorthEast');
    else
        hleg1=legend(num2str(T_double(1)),num2str(T_double(2)),num2str(T_double(3)),'Location','NorthEast');
    end
    hold off;
end

for no_plt = 3:1:7
    switch no_plt
        case 3
            axes(handles.axes3);
        case 4
            axes(handles.axes4);
        case 5
            axes(handles.axes5);
        case 6
            axes(handles.axes6);
        case 7
            axes(handles.axes7);
    end
    cla  
    hold all;
    for i=1:1:ii
    plot(tm, yy(:,no_plt,i),lineopt(i).linetype,'LineWidth',2);
        if on_semilogY(no_plt) == 1
        set(gca,'yscale','log') %use for semilogy
        else
          set(gca,'yscale','linear')  
        end
    end
    xlabel('Time [days]');ylabel(plotparam(no_plt).ytitle);
    hleg1=legend(num2str(T_double(1)),num2str(T_double(2)),num2str(T_double(3)),'Location','NorthEast');
    hold off;
end

%save the data
fid1 = fopen('calcdata.txt','w');
fprintf(fid1,'Time [days], Total volumes [cm^3] for three T_double\n');
for i =1:1:size(T)
fprintf(fid1,'%5.2f %8.4e %8.4e %8.4e\n', tm(i), yy(i, 1, 1), yy(i, 1, 2), yy(i, 1, 3));
end

% Tumor volume doubling time at GKSRS
tgk = rad_on(1);
for i=1:ii
td_gksrs(i) = yy(tgk,6,i)/2;
end
fprintf(fid1,'\nDoubling time at GKSRS [days]\n');
fprintf(fid1,'%8.4f %8.4f %8.4f\n', td_gksrs(1), td_gksrs(2), td_gksrs(3));

if (cont_kill_model == 1)
fprintf(fid1,'\n%25s %12.4f\n', 'Cell kill probability = ', q);
fprintf(fid1,'The total volume of tumor cells killed by radiation [volume ratio/integ of kill term]\n');
fprintf(fid1,'%12.6f %12.6f %12.6f\n', total_cellkill(1), total_cellkill(2), total_cellkill(3));
fprintf(fid1,'%12.6f %12.6f %12.6f\n', total_cellkill2(1), total_cellkill2(2), total_cellkill2(3));
fprintf(fid1,'The cell survival fraction of tumor cells after irradiation\n');
fprintf(fid1,'%12.3e %12.3e %12.3e\n', cellsurvival_rate(1), cellsurvival_rate(2), cellsurvival_rate(3));
end
fprintf(fid1,'\n%25s %12.3e','LQ_survival fraction = ',exp(-chi));

fclose(fid1);

% --------------------------------------------------------------------
function tagfile_Callback(hObject, eventdata, handles)
% hObject    handle to tagfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function tag_snapshot_Callback(hObject, eventdata, handles)
% hObject    handle to tag_snapshot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = getframe(handles.figure1);
imwrite(I.cdata,'snapshot.jpg');

% --------------------------------------------------------------------
function tag_exit_Callback(hObject, eventdata, handles)
% hObject    handle to tag_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figure1);



function eta_vasc_Callback(hObject, eventdata, handles)
% hObject    handle to eta_vasc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eta_vasc as text
%        str2double(get(hObject,'String')) returns contents of eta_vasc as a double


% --- Executes during object creation, after setting all properties.
function eta_vasc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eta_vasc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_repair_Callback(hObject, eventdata, handles)
% hObject    handle to T_repair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_repair as text
%        str2double(get(hObject,'String')) returns contents of T_repair as a double


% --- Executes during object creation, after setting all properties.
function T_repair_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_repair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Y_initial_Callback(hObject, eventdata, handles)
% hObject    handle to Y_initial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Y_initial as text
%        str2double(get(hObject,'String')) returns contents of Y_initial as a double


% --- Executes during object creation, after setting all properties.
function Y_initial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Y_initial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function tag_help_Callback(hObject, eventdata, handles)
% hObject    handle to tag_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function about_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('TUMORSIM Version 5.5: Last update on 12/13/2015');


function on_semilogy_Callback(hObject, eventdata, handles)
% hObject    handle to on_semilogy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of on_semilogy as text
%        str2double(get(hObject,'String')) returns contents of on_semilogy as a double


% --- Executes during object creation, after setting all properties.
function on_semilogy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to on_semilogy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function tagDataIn_Callback(hObject, eventdata, handles)
% hObject    handle to tagDataIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%import a set of experimental data of tumor growth from .txt file
%1st column = time in day
%2nd column = volume in cm^3
global dose p q rad_on expdatain voldata eta_clearance
fid = fopen('volume_data.txt');
tempin = fscanf(fid,'%f %f');
fclose(fid);
nsets = size(tempin,1)/2;
voldata = zeros(nsets,2);
for i=1:1:nsets
    voldata(i,1)=tempin((i-1)*2+1);
    voldata(i,2) = tempin(i*2);
end
expdatain = 1;


% --- Executes on button press in pltexpdata.
function pltexpdata_Callback(hObject, eventdata, handles)
% hObject    handle to pltexpdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pltexpdata


% --- Executes on button press in Optparams.
function Optparams_Callback(hObject, eventdata, handles)
% hObject    handle to Optparams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%estimate unknown model parameters in ODEs
%used with TUMORSIM
%written by Yoichi Watanabe
%based on odeParamID of MathWork
%last revised on 11/19/2013
%c(1) = lammda
%c(2) = mu
%c(3) = theta
%y(1) = tumor volume
%y(2) = volume of non-dividing cells
%y(3) = tumor growth rate
global p q rad_on expdatain voldata eta_clearance instant_kill_model
%=================== input parameters ======
dose = str2num(get(handles.dose,'String')); %single dose [Gy]
dt = str2num(get(handles.T_step,'String')); %time step for simulation [day]
T_duration = str2num(get(handles.T_duration,'String')); %the total time in days
tmax = T_duration(2);
Y0 = str2num(get(handles.Y_initial,'String')); %initial condition
rad_effective = str2num(get(handles.rad_on,'String')); %hours when the radiation is killing the cells
T_star = str2num(get(handles.T_cycle,'String')); %day cell cycle time
c00 = str2num(get(handles.paramguess,'String'));
rad_on = rad_effective*(1/T_star); %turn on radiation, single instant dose
vol_norm = 1.0;
alog2 = log(2);
%% load "experimental data"
fid = fopen('volume_data.txt');
tempin = fscanf(fid,'%f %f');
fclose(fid);
nsets = size(tempin,1)/2;
voldata = zeros(nsets,2);
for i=1:1:nsets
    expTime(i)=tempin((i-1)*2+1); %non-dimensional time
    expY(i) = tempin(i*2); %tumor volume in cm^3 divided by one cm^3
end

figure('Name','Parameter Optimization','NumberTitle','off') 
semilogy(expTime, expY, 'or')
%non-dimensionalize Time and Y
expTime0 = expTime/T_star;
expY0 = expY/vol_norm;
%% ODE Information
%initial guess of paramters
c0(1) = alog2/c00(1)*T_star;
c0(2) = alog2/c00(2)*T_star;
c0(3) = c00(3);
Y0(3) = c0(1);
%----------------------------
i1 = floor(T_duration(2)/dt);
for j=1:1:i1
    tt(j) = (dt/T_star)*(j-1);
end
tSpan(1) = 0.0;
tSpan(2) = tt(i1);
%initial condition for ODEs
%% Initial guess
[T, Y] = ode15s(@(t,y)fun_vdp9a(t,y,c0),tSpan,Y0);
%YY = Y(:,1)+Y(:,2);
YY = Y(:,1);
TT = T*T_star;
hold on
semilogy(TT, YY, '-g')

%% Set up optimization
myObjective = @(x) objFcn(x, expTime, expY, tt, Y0);
lb = [0 0 0];
ub = [10 10 10];

bestc = lsqnonlin(myObjective, c0, lb, ub);
%% Plot best result
[T, Y] = ode15s(@(t,y)fun_vdp9a(t,y,bestc),tSpan,Y0);
%YY = Y(:,1)+Y(:,2);
YY = Y(:,1);
TT = T*T_star;
hold on
semilogy(TT, YY, '-b')
legend('Exp Data','Initial Param','Best Param');
xlabel('Time [days]');
ylabel('Tumor volume [cm^3]');

Td_opt = alog2/bestc(1)*T_star;
Tr_opt = alog2/bestc(2)*T_star;
Theta = bestc(3);
optresults = sprintf('%7.2f %7.2f %7.4f',Td_opt, Tr_opt, Theta);
set(handles.optresults,'String',optresults);

function paramguess_Callback(hObject, eventdata, handles)
% hObject    handle to paramguess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of paramguess as text
%        str2double(get(hObject,'String')) returns contents of paramguess as a double


% --- Executes during object creation, after setting all properties.
function paramguess_CreateFcn(hObject, eventdata, handles)
% hObject    handle to paramguess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function optresults_Callback(hObject, eventdata, handles)
% hObject    handle to optresults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of optresults as text
%        str2double(get(hObject,'String')) returns contents of optresults as a double


% --- Executes during object creation, after setting all properties.
function optresults_CreateFcn(hObject, eventdata, handles)
% hObject    handle to optresults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%----------------------------------------------------------------------------
function dy = fun_vdp9a(t, y, c)
global p q rad_on expdatain voldata eta_clearance instant_kill_model

if instant_kill_model == 0
    if (t >= rad_on(1) && t <= rad_on(2))
        pp = p;
        qq = q;
        cc1 = 0.0;
    else
        pp = 1.0;
        qq = 0.;
       cc1 = c(1);
%        cc1 = 1.0;
    end
else
        pp  = 1.0;
        qq = 0.0;
        cc1 = c(1);
end
dy= zeros(3,1);
dy(1)=2*y(3)*pp*y(1) - qq*y(1) + c(2)*y(2);
dy(2)=qq*y(1) - (c(2)+eta_clearance)*y(2);
dy(3)=(c(3)-1.0)*cc1*y(3);
%dy(3)=(c(3)-1.0)*cc1*y(3)^2;

function cost = objFcn(x,expTime,expY,tSpan,z0)
ODE_Sol = ode15s(@(t,z)fun_vdp9a(t,z,x), tSpan, z0);
simY = deval(ODE_Sol, expTime);
%simYY = simY(1,:) + simY(2,:);
simYY = simY(1,:);
cost = simYY - expY;

function T_clear_Callback(hObject, eventdata, handles)
% hObject    handle to T_clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_clear as text
%        str2double(get(hObject,'String')) returns contents of T_clear as a double


% --- Executes during object creation, after setting all properties.
function T_clear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kappa_Callback(hObject, eventdata, handles)
% hObject    handle to kappa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kappa as text
%        str2double(get(hObject,'String')) returns contents of kappa as a double


% --- Executes during object creation, after setting all properties.
function kappa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kappa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in instantkill.
function instantkill_Callback(hObject, eventdata, handles)
% hObject    handle to instantkill (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of instantkill


% --- Executes on button press in cont_cell_kill.
function cont_cell_kill_Callback(hObject, eventdata, handles)
% hObject    handle to cont_cell_kill (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cont_cell_kill

% ---- calculate the total number of cells killed by radiation for Model B
function total_cellkill = intgl_killterm(rad_on, T, y, p, dt)
jj = size(T, 1);
cc = 0.0;
for i = 1:1:jj
    t = T(i);
    if (t >= rad_on(1) & t <= rad_on(2))
        cc = cc + p*dt*y(i);
    end
end
total_cellkill = cc;
