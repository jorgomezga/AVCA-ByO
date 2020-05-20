function [x_intersect,y_intersect]=intersectcurves(X, Y1, Y2)
% Intersection point for cumulative curves
% 
% INPUT:
%       X - x axis values
%       Y1, Y2 - cumulative curves
%
% OUTPUT:
%       x_intersect, y_intersect: intersection point (or its interpolation)

vDiff=Y2-Y1;
if vDiff(1)<0
    vDiff=-vDiff;
end

ind=find(vDiff==0);

if ~isempty(ind)
    x_intersect=X(ind);
    y_intersect=Y1(ind);
else
    ind=find(vDiff<0);
    ind=ind(1);
    x_intersect=X(ind-1)+(X(ind)-X(ind-1))/2;
    y_intersect=Y1(ind-1)+abs((Y1(ind)-Y1(ind-1))/2);
    
end

