function H = BarPlot(X)
% =======================================================================
% Creates a bar graph with positive data values stacked on the positive
% quadrant and negative data values stacked on the negative quadrant
% =======================================================================
% H = BarPlot(X)
% -----------------------------------------------------------------------
% INPUT
%   - X: data to plot [nobs x nvars]
% -----------------------------------------------------------------------
% OUTPUT
%   - H: handle to graph
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com
colormap winter
H(1,:) = bar((X).*(X>0),'stacked'); 
for ii=1:size(H,2)
    H(1,ii).EdgeColor = H(1,ii).FaceColor;
end
hold on;
H(2,:) = bar((X).*(X<0),'stacked');
for ii=1:size(H,2)
    H(2,ii).EdgeColor = H(2,ii).FaceColor;
end
