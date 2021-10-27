function anaKine
% Controls the acquisition of kinematics for the lionfiosh project
% NOTE: Defaults to batch mode, if no sequence information provided


%% Process inputs

% Define data paths
datapath = 'C:\Users\anpet\Documents\Repositories\Teaching\exData\S18';%uigetdir;

% Define path to save analysis
selpath = 'C:\Users\anpet\Documents\Repositories\Teaching\exOutput';%uigetdir;

% Load raw datafiles
load([datapath filesep 'Body.mat']);
load([datapath filesep 'initial conditions.mat']);
load([datapath filesep 'Landmarks.mat']);

% Initilize
ang = 0;

%Parameters
sp = [];
d  = [];



%% Nonspecific Variables

% Start frame (if need to start from nonzero frame)
d.startFrame   = 2;
d.strikeFrame  = Body.strike_frame;

% vars for trimming frames
c = d.startFrame;
e = find(Body.frames == Body.strike_frame);


% if necessary, correct strike frame
if isempty(e)
    warning('corrected strike frame is beyond analyzed frame: need to fix')
    save([datapath filesep 'badendFrame.mat'],'e')
    e = find(Body.frames == Body.frames(end));;
    d.strikeFrame = Body.frames(end);
end

% conversion factor (pixels/m)
d.cF = iC.calconst;

% Video height to calibrate Y data
V.height = iC.vHeight;  % v.Height

% Time
d.t = (Body.frames-min(Body.frames))./iC.frameRate;
d.t = d.t(c:e);

% Initial period used to identify initial heading (s)
d.startPeriod = 0.5;

% Store away general data
d.frame        = Body.frames(c:e);
d.framerate    = iC.frameRate;
d.strike_frame = length(Body.frames(c:e));



%% Initilize and Adjust Fish Variables


% Positional

% Predator trajectories (in pixels)
d.xPd_pix = Body.x(c:e,1);
d.yPd_pix = V.height - Body.y(c:e,1);

% Prey trajectories (in pixels)
d.xPy_pix = Body.x(c:e,2);
d.yPy_pix = 

% (in meters)
% Predator
d.xPd = 
d.yPd = 

% Prey
d.xPy = 
d.yPy = 



%% Analyze coordinate data

% Start angle
d.startAng = calcdata('start angle',d);

% Extract spline fits from data
d.sp = calcdata('splines',d);

% Calculate speeds
[d.sp.spdPy,d.sp.spdPd,d.sp.dist] = calcdata('speed',d.sp);


% Calculate angles
[d.sp.bearing,d.sp.alpha,d.sp.thetaPy,d.sp.thetaPd] = ...
    calcdata('angles',d);



%% Calculate data

% Distance (hypot)



%% Plot




%% Save Output

% Save data file
save([selpath filesep 'anaSeq.mat'],'d');

display('Analysis complete')




function varargout = calcdata(action,dIn)
%% Calculate data


if strcmp(action,'angles') % Angles ---------------------------------------
    
    % Run calculation
    [bearing,alpha,thetaPy,thetaPd] = calc_angles(dIn);

    
    % Store output
    varargout{1} = bearing;
    varargout{2} = alpha;
    varargout{3} = thetaPy;
    varargout{4} = thetaPd;

elseif strcmp(action,'start angle') % Start angle -------------------------
    % Fit to initial points
    cXprey = polyfit(dIn.t,dIn.xPy,1);
    cYprey = polyfit(dIn.t,dIn.yPy,1);
    
    % Visualize
    if 0
        figure;
        plot(dIn.xPy,dIn.yPy,'r.',...
            polyval(cXprey,dIn.t),polyval(cYprey,dIn.t),'k-')
        axis equal
    end
    
    % Initial heading of prey
    varargout{1} = atan2(polyval(cYprey,dIn.startPeriod)-polyval(cYprey,0),...
        polyval(cXprey,dIn.startPeriod)-polyval(cXprey,0));
    
elseif strcmp(action,'speed') % Speed ------------------------------------
    
    % Run calculation
    [spdPy,spdPd,dist] = calc_spd(dIn);
    
    % Outputs
    varargout{1} = spdPy;
    varargout{2} = spdPd;
    varargout{3} = dist;
    
elseif strcmp(action,'splines') % Splines ---------------------------------
    
    varargout{1} = calc_splines(dIn);
    
else
    error('action not identified')
end

function sp = calc_splines(d)
%% Calculate splines

% Whether to plot data
spCheck = 0;

frame  = d.frame';
xPrey  = d.xPy;
yPrey  = d.yPy;
xPred  = d.xPd;
yPred  = d.yPd;

% generate a time vector
sp.t = d.t;

% Store frame numbers
sp.frStart = d.startFrame;
sp.frEnd   = length(frame);


% Smooth prey data
xPrey_sm = smooth(xPrey);
yPrey_sm = smooth(yPrey);

% Spline fit prey position
sp.xPy = spline(sp.t,xPrey_sm);
sp.yPy = spline(sp.t,yPrey_sm);

sp.xPrey_fit = fnval(sp.t,sp.xPy);
sp.yPrey_fit = fnval(sp.t,sp.yPy);

% Smooth predatory data data
xPred_sm = smooth(xPred);
yPred_sm = smooth(yPred);

% Spline fit predator rostrum position
sp.xPd = spline(sp.t,xPred_sm);
sp.yPd = spline(sp.t,yPred_sm);

sp.xPred_fit = fnval(sp.t,sp.xPd);
sp.yPred_fit = fnval(sp.t,sp.yPd);

function [bearing,alpha,thetaPy,thetaPd] = calc_angles(d)
%% Calculate angular quantities from trajectory data

% Solve for discrete points from splines
xPrey_fit = fnval(d.sp.t,d.sp.xPy);
yPrey_fit = fnval(d.sp.t,d.sp.yPy);
xPred_fit = fnval(d.sp.t,d.sp.xPd);
yPred_fit = fnval(d.sp.t,d.sp.yPd);

% % Length of data vectors
% numData = length(xPred_fit);

% Vector from predator to prey 'baseline vector'
vBase = [xPrey_fit - xPred_fit,...
    yPrey_fit - yPred_fit];


% Bearing angle using pred velocity for heading angle ---------------------

% Prey velocity from spline
spdxPy = fnder(d.sp.xPy);
spdyPy = fnder(d.sp.yPy);
dxPrey_fit = fnval(d.sp.t,spdxPy);
dyPrey_fit = fnval(d.sp.t,spdyPy);

% Predator velocity from spline
d.sp.spdxPd = fnder(d.sp.xPd);
d.sp.spdyPd = fnder(d.sp.yPd);
dxPred_fit = fnval(d.sp.t,d.sp.spdxPd);
dyPred_fit = fnval(d.sp.t,d.sp.spdyPd);

dxPred_fitsm = smooth(dxPred_fit);
dyPred_fitsm = smooth(dyPred_fit);
dxPrey_fitsm = smooth(dxPrey_fit);
dyPrey_fitsm = smooth(dyPrey_fit);

d.sp.dyPred_fit = spline(d.sp.t,dyPred_fitsm);
d.sp.dxPred_fit = spline(d.sp.t,dxPred_fitsm);
d.sp.dyPrey_fit = spline(d.sp.t,dxPrey_fitsm);
d.sp.dxPrey_fit = spline(d.sp.t,dyPrey_fitsm);

dyPred_fitsm = fnval(d.sp.t,d.sp.dyPred_fit);
dxPred_fitsm = fnval(d.sp.t,d.sp.dxPred_fit);
dyPrey_fitsm = fnval(d.sp.t,d.sp.dyPrey_fit);
dxPrey_fitsm = fnval(d.sp.t,d.sp.dxPrey_fit);

% Predator heading angle using velocity
thetaPd = atan2(dyPred_fitsm,dxPred_fitsm);


% Angle of baseline vector (alpha)
alpha = atan2(vBase(:,2),vBase(:,1));

% Bearing angle
bearing = atan2(sin(alpha - thetaPd),cos(alpha - thetaPd));



% Speed of prey
vPrey = sqrt(sum([dxPrey_fitsm, dyPrey_fitsm].^2,2));

% Speed of predator
vPred = sqrt(sum([dxPred_fitsm, dyPred_fitsm].^2,2));

% Prey heading angle using velocity
thetaPy= atan2(dyPrey_fit,dxPrey_fit);



% Angle between prey velocity and baseline vector
beta= atan2(sin(alpha - thetaPy),cos(alpha - thetaPy));


function [spdPy,spdPd,dist] = calc_spd(sp)
%% Return speeds for predator and prey, based on spline data

% Coordinates
xPy  = fnval(sp.xPy,sp.t);
yPy  = fnval(sp.yPy,sp.t);
xPd  = fnval(sp.xPd,sp.t);
yPd  = fnval(sp.yPd,sp.t);

% Distance between
dist = hypot(xPy-xPd,yPy-yPd);

% Derivatives along x and y
dXPy  = fnval(fnder(sp.xPy),sp.t);
dYPy  = fnval(fnder(sp.yPy),sp.t);

% Speed
spdPy = hypot(dXPy,dYPy);

% Derivatives along x and y
dXPd  = fnval(fnder(sp.xPd),sp.t);
dYPd  = fnval(fnder(sp.yPd),sp.t);

% Speed
spdPd = hypot(dXPd,dYPd);





