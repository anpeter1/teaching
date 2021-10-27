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

% % Frame of reference angle (in Degrees)
% d.FoRang = L.ang_pd;
% d.FoRfrm = L.FoRframes;

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
d.yPy_pix = V.height - Body.y(c:e,2);

% (in meters)
% Predator
d.xPd = d.xPd_pix*d.cF;
d.yPd = d.yPd_pix*d.cF;

% Prey
d.xPy = d.xPy_pix*d.cF;
d.yPy = d.yPy_pix*d.cF;



% % Predator body orientation from tform
% 
% % extract tform
% for f = c:length(Body.Rotation.tform)
%     tform(f) = Body.Rotation.tform(f);
%     d.tform.Pd.x(f,1) = tform(f).T(1,1);
%     d.tform.Pd.y(f,1) = tform(f).T(1,2);
% end
% clear f
% 
% % correct for c start
% d.tformPd.x(:,1) = d.tform.Pd.x(c:e,1);
% d.tformPd.y(:,1) = d.tform.Pd.y(c:e,1) ;



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

% [d.sp.bearing,d.sp.alpha,d.sp.thetaPy,d.sp.thetaPd,d.sp.thetaPdH,d.sp.bearingH] = ...
%     calcdata('angles',d);



%% Calculate data

% Distance
d.dist = hypot(d.xPy-d.xPd,d.yPy-d.yPd);
sp.dist = hypot(d.sp.xPrey_fit-d.sp.xPred_fit, d.sp.yPrey_fit-d.sp.yPred_fit);


%% Plot

if 0
    
    if plotsplineeval
        figure
        title([num2str(seqNum),'   Fish ID'])
        subplot(2,1,1)
        hold on
        plot(D.t,D.xPy,'bo',D.t,D.yPy,'bo')
        fnplt(sp.xPrey,'y')
        fnplt(sp.yPrey,'y')
        title('Prey Splines')
        %         xlabel('x')
        %         ylabel('Time')
        subplot(2,1,2)
        hold on
        plot(D.t,D.xPd,'ro',D.t,D.yPd,'ro')
        fnplt(sp.xPred,'b')
        fnplt(sp.yPred,'b')
        title('Pred Splines')
        %         xlabel('x')
        %         ylabel('Time')
        
        
        figure
        
        subplot(2,1,1)
        hold on
        plot(D.t,d.spdPy,'bo')
        fnplt(sp.spdPy,'g')
        title(['Prey Splines ',num2str(sp.tol.py),' tol,'])
        subplot(2,1,2)
        hold on
        plot(D.t,d.spdPd,'ro')
        fnplt(sp.spdPd,'g')
        title(['Pred Splines ',num2str(sp.tol.pd),' tol,'])
        
    end
    
    if plotposition
        figure
        hold on
        plot(d.xTank,d.yTank)
        plot(sp.xPyVal,sp.yPyVal,'b')
        plot(sp.xPdVal,sp.yPdVal,'r')
        set(gca, 'YDir','reverse')
        title('Fish Trajectories')
        
    end
    
    if plotdist
        figure
        hold on
        plot(D.t,d.dist,'gx')
        plot(D.t,sp.dist,'k')
        title('Distance btw Pd/Py')
        xlabel('Time (s)')
        ylabel('Distance (m)')
    end
    
    if plotspeed
        figure
        subplot(2,1,1)
        hold on
        plot(D.t(bump:end),d.spdPd(bump:end),'bo',D.t(bump:end),d.spdPy(bump:end),'bo')
        fnplt(sp.xPrey,'y')
        fnplt(sp.yPrey,'y')
        title('Prey Splines')
        %         xlabel('x')
        %         ylabel('Time')
        subplot(2,1,2)
        hold on
        plot(D.t,D.xPd,'ro',D.t,D.yPd,'ro')
        fnplt(sp.xPred,'b')
        fnplt(sp.yPred,'b')
        title('Pred Splines')
        
        
        figure
        %velocity
        set(0,'DefaultFigureWindowStyle','docked')
        subplot(2,1,1)
        hold on
        plot(D.t(bump:end),sp.spdPyVal(bump:end),'b')
        plot(D.t(bump:end),sp.spdPdVal(bump:end),'r')
        %         ylim([0 0.2])
        %         ylim([0 2])
        title('Fish Velocity')
        xlabel('Time (s)')
        ylabel('Velocity (m/s)')
        
        
        %Velocity Raw
        subplot(2,1,2)
        hold on
        plot(D.t(bump:end),d.spdPd(bump:end),'r')
        plot(D.t(bump:end),d.spdPy(bump:end),'b')
        %         ylim([0 1])
        %         ylim([0 2])
        title('Raw Velocity')
        xlabel('Time (s)')
        ylabel('Velocity (m/s^2)')
        
        
        % Accel plots
        if plotacc
            figure
            %acceleration
            subplot(2,1,2)
            hold on
            plot(D.t,sp.accPdVal,'r')
            plot(D.t,sp.accPyVal,'b')
            ylim([0 1])
            title('Fish Acceleration')
            xlabel('Time (s)')
            ylabel('Velocity (m/s^2)')
        end
        %
        %     figure
        %     plot(D.t,d.spdPy,'b',D.t,d.spdPd,'r')
        %         title('Fish velocity')
        %         xlabel('Time (s)')
        %         ylabel('Velocity (m/s)')
        %
    end
    
    % Angle plots
    if plotang
        figure
        hold on
        plot(D.t,D.angPy,'bo',D.t,D.angPd,'ro')
        plot(D.t,sp.angPyVal,'g',D.t,sp.angPdVal,'k')
        
        figure
        subplot(2,2,1)
        plot(D.t,d.thetaPy,'b',D.t,d.thetaPd,'r')
        
        subplot(2,2,2)
        plot(D.t,d.bearing)
        
        subplot(2,2,3)
        plot(D.t,d.alpha)
        
        subplot(2,2,4)
        plot(D.t,d.beta)
        
    end
    
    
end


%% Save Output

% Save data file
save([selpath filesep 'anaSeq.mat'],'d');

display('Analysis complete')




function varargout = calcdata(action,dIn)
%% Calculate data


if strcmp(action,'angles') % Angles ---------------------------------------
    
    % Run calculation
    [bearing,alpha,thetaPy,thetaPd] = calc_angles(dIn);
%     [bearing,alpha,thetaPy,thetaPd,thetaPdH,bearingH] = calc_angles(dIn);
    
    % Store output
    varargout{1} = bearing;
    varargout{2} = alpha;
    varargout{3} = thetaPy;
    varargout{4} = thetaPd;
%     varargout{5} = thetaPdH;
%     varargout{6} = bearingH;
%     
    
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
% function [bearing,alpha,thetaPy,thetaPd,thetaPdH,bearingH] = calc_angles(d)
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

% % Predator heading from tform
% thPdH = atan2(d.tformPd.y,d.tformPd.x);
% 
% % Calibrate theta from heading with FoR angles
% 
% for f = 1:length(d.FoRfrm)
%     aRot(f,1) = thPdH(find(d.frame==d.FoRfrm(f)));
% end
% clear f
% 
% aMeas = d.FoRang;
% aDiff = unwrap(aRot) + unwrap(aMeas);
% aFinal = thPdH - mean(aDiff);
% thetaPdH = aFinal;



% diff = d.FoRang*pi/180;
% thetaPdH = -diff + thPdH(:);
% thetaPdH = thetaPdH*180/pi;

% Angle of baseline vector (alpha)
alpha = atan2(vBase(:,2),vBase(:,1));

% Bearing angle
bearing = atan2(sin(alpha - thetaPd),cos(alpha - thetaPd));

% % Bearing angle from tform theta
% bearingH = atan2(sin(alpha - thetaPdH),cos(alpha - thetaPdH));

% Speed of prey
vPrey = sqrt(sum([dxPrey_fitsm, dyPrey_fitsm].^2,2));

% Speed of predator
vPred = sqrt(sum([dxPred_fitsm, dyPred_fitsm].^2,2));

% Prey heading angle using velocity
thetaPy= atan2(dyPrey_fit,dxPrey_fit);

% % Prey heading angle from tform
% thPyH = atan2(d.tformPy.y,d.tformPy.x);

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





