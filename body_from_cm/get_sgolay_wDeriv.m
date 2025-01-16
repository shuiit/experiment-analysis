function [SG0, SG1, SG2] = get_sgolay_wDeriv(y, N, F, FPS)
% FPS = 1/dt
%{
sgDegree = 5 ;
sgWindow = 21 ;
halfWin  = floor(sgWindow/2) ;

[xfilt, vx, ~] = get_sgolay_wDeriv(x_raw, sgDegree, sgWindow, FPS);
[yfilt, vy, ~] = get_sgolay_wDeriv(y_raw, sgDegree, sgWindow, FPS);
%}
halfWin  = ((F+1)/2) - 1;
[~,g] = sgolay(N,F);

sz = length(y) ;
SG0 = zeros(sz,1) ;
SG1 = zeros(sz,1) ;
SG2 = zeros(sz,1) ;

for n = (F+1)/2 : sz-(F+1)/2
    % Zeroth derivative (smoothing only)
    SG0(n) = dot(g(:,1),y(n - halfWin:n + halfWin));
    % 1st differential
    SG1(n) = dot(g(:,2),y(n - halfWin:n + halfWin));
    % 2nd differential
    SG2(n) = 2*dot(g(:,3)',y(n - halfWin:n + halfWin))';
end

SG1 = SG1 * FPS ;
SG2 = SG2 * (FPS^2) ;

end