function [f, fLift, fDrag] = quasiSteadtForceSimple_mk2(vbs, geom)
% vbs is a 6-componenet velocity vector of the wing velocity in the WING
% frame of reference
% vbs(1:3) = angular velocity in rad/sec [omega_x ; omega_y ; omega_z] ;
% vbs(4:6) = translational velocity of the entire fly body.

% find tip velocity in the wing frame, by combining rotational and
% translational velocities of the wing

%vbs(1:3) = 0 ;

rv    = zeros(4, 1) ;
rv(1) = geom.span / 2; % previously was span/2.
rv(4) = 1 ;

v = wedge(vbs) * rv ;  % v = cross(omega, span) + translational velocity

vyp = v(2) ; % velocity component perpendicular to wing surface
vzp = v(3) ; % velocity component parallel to chord vector

U     = [ 0, vyp, vzp] ;
U2    = (vyp.^2 + vzp.^2) ;

if (U2==0)
    %if zero velocity return zero force
    f = zeros(3,1) ;
    fLift = zeros(3,1) ;
    fDrag = zeros(3,1) ;
    return ;
end
disp(sqrt(U2)) ;

Uhat  = U / sqrt(U2) ;
alpha = atan2(vzp, vyp) ; 

airDensity = geom.airDensity ; 

r22 = geom.r22 ; % 0.4 ; % 0.3 ; % geometric factor, second moment of area
wingArea = geom.span * geom.chord * pi / 4; % area of an ellipse

CLMAX = 1.8 ; % see Dickinson Science 1999 and Whitney & Wood JFM 2010
CDMAX = 3.4 ;
CD0   = 0.4 ;
%CROT  = 1.55 ; %#ok<NASGU>
%rotShapeIntegral = 0.5 ; %#ok<NASGU>

% lift and drag coefficients
%CD = (CDMAX/2) * ( 1 - cosd(2*alpha)) ;
CD = (CDMAX+CD0)/2 - ((CDMAX-CD0)/2) * cos(2*alpha) ;
CL = CLMAX * abs(sin(2*alpha)) ;

fDrag = -0.5 * CD * airDensity * wingArea * r22 * Uhat * U2 ;

%chordHat = [0 1 0 ] ;
%altSpanVec = cross(Uhat, chordHat) ;
%altSpanVec = altSpanVec / norm(altSpanVec) ;
%Lhat = cross(Uhat, altSpanVec) ;
spanHat = [1 0 0 ] ; 
Lhat = cross(Uhat, spanHat) ; % perpendicular to Uhat
Lhat = Lhat / norm(Lhat) ; % but still need to determine sign

% Lhat and Uhat should be on different sides on the wing, so they should
% have different signs of their zprime component

% sign(Lhat(3)) * sign(Uhat(3)) is negative if Lhat is on the right
% direction. Hence we have:

Lsign = -sign(Lhat(3)) * sign(Uhat(3)) ;

Lhat  = Lhat * Lsign ;

fLift = + 0.5 * CL * airDensity * wingArea * r22 * Lhat * U2 ;

%if (fLift(3)<0) ; keyboard ; end ;

f = fDrag + fLift ;
f=f' ;

%{
if (includeRotationalForce)
    
    % rotational force
    % ----------------
    % see eq. 11 in Sane and Dickinson JEB 2001
    % calculate magnitude and direction separately
    
    currt = t(it) ;
    etadot_rad_sec = eval_fd(currt, fdEtaR_lab, 1) * (pi/180) * 1000 ; % 1000 from msec to sec
    % NEED TO USE OMEGA=psidot + phidot*sin(theta)
    
    psiDot = 
    
    spanHat = data.rightSpanHats(it,:);
    rotHat = cross(spanHat, chordHat) * sign(etadot_rad_sec) ;
    rotHat = rotHat/norm(rotHat) ;
    
    frot = CROT * airDensity * sqrt(U2) * abs(etadot_rad_sec) * ...
        wingChordMeters^2 * wingSpanMeters * rotShapeIntegral ;
    
    FrotR(it,:) = frot * rotHat ;
    
    
end

%}
return



