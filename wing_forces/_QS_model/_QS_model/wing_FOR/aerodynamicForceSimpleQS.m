function [f, fLift, fDrag] = aerodynamicForceSimpleQS(u, geom)
% taken from dpsi_mk3.m (folder springy)
% 
% inputs
% ------
% u is the wing tip velocity in the wing (!!) frame of reference (3-vector)
%
% wing frame of reference is: 
%   x == span
%   y == chord
%   z == cross(span, chord)
% make sure you handle left and right wing correctly...
%
% geom is a struct with the following 4 fields:
%
% geom.airDensity
% geom.r22          % 0.4 ; % 0.3 ; % geometric factor, second moment of area
% geom.span
% geom.chord
%
% outputs
% -------
% f     -- total force in the wing frame of reference (vector)
% fLift -- lift force (vector vertical to u)
% fDrag -- drag force (vector parallel to u)


% find angle of attack

uy = u(2) ;
uz = u(3) ;
% no ux since quasi-steady does not consider ux

U     = [ 0 uy uz] ;
U2    = (uy.^2 + uz.^2) ;

if (U2==0)
    %if zero velocity return zero force
    f = zeros(3,1) ;
    fLift = zeros(3,1) ;
    fDrag = zeros(3,1) ;
    return ;
end

Uhat  = U / sqrt(U2) ;
alpha = atan2(uz, uy) ; 

airDensity = geom.airDensity ; 

r22 = geom.r22 ; % 0.4 ; % 0.3 ; % geometric factor, second moment of area
wingArea = geom.span * geom.chord * pi / 4; % area of an ellipse

CLMAX = 1.8 ; % see Dickinson Science 1999 and Whitney & Wood JFM 2010
CDMAX = 3.4 ;
CD0   = 0.4 ;
%CROT  = 1.55 ;
%rotShapeIntegral = 0.5 ; 

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
% have different signs of their zp rime component

% sign(Lhat(3)) * sign(Uhat(3)) is negative if Lhat is on the right
% direction. Hence we have:

Lsign = -sign(Lhat(3)) * sign(Uhat(3)) ;

Lhat  = Lhat * Lsign ;

fLift = + 0.5 * CL * airDensity * wingArea * r22 * Lhat * U2 ;

f = fDrag + fLift ;
f=f' ;

end
