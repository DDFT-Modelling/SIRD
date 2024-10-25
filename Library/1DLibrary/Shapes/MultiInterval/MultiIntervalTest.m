geom1.N = 20;
geom1.yMin = -2;
geom1.yMax = 0;

geom2.N = 10;
geom2.yMin = 0;
geom2.yMax = 1;

intervals(1).interval = 'SpectralLine';
intervals(2).interval = 'SpectralLine';
intervals(1).geom = geom1;
intervals(2).geom = geom2;

BCs = struct;
BCs(1,2).function = 'BCmatch';

MI = MultiInterval(intervals);
MI.ComputeUniformPts;


Pts = MI.Pts;
Ind = MI.Ind;
Diff = MI.Diff;
Int = MI.Int;
nIntervals = MI.nIntervals;
Intersections = MI.Intersections;

y = Pts.y;

Dy = Diff.Dy;

% MI.Plot(y);
% MI.Plot(y.^2);

max(abs(Dy*(y) - 1))
max(abs(Dy*(y.^2) - 2*y))

max(abs(Int*y - (MI.Pts.yMax^2 - MI.Pts.yMin^2)/2 ))
max(abs(Int*(y.^2) - (MI.Pts.yMax^3 - MI.Pts.yMin^3)/3))

% figure
% MI.PlotBound({'left','right'});