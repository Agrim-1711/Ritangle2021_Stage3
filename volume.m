clc;
%define the 187 spheres in a matrix array [x y z R]
Points= [8 8 8 8
15 8 8 1
7 8 8 7
7 14 8 1
7 7 8 6
12 7 8 1
6 7 8 5
6 11 8 1
6 6 8 4
9 6 8 1
5 6 8 3
5 8 8 1
5 5 8 2
4 5 8 1
6 5 8 1
2 2 2 2
14 2 2 2
2 14 2 2
14 14 2 2
2 2 1 1
14 2 1 1
2 14 1 1
14 14 1 1
2 7 1 1
14 7 1 1
2 9 1 1
14 9 1 1
7 2 1 1
9 2 1 1
7 14 1 1
9 14 1 1
5 1 3 1
11 1 3 1
5 15 3 1
11 15 3 1
1 5 3 1
15 5 3 1
1 11 3 1
15 11 3 1
1 3 5 1
15 3 5 1
1 13 5 1
15 13 5 1
3 1 5 1
13 1 5 1
3 15 5 1
13 15 5 1
3 5 1 1
13 5 1 1
3 11 1 1
13 11 1 1
5 3 1 1
11 3 1 1
5 13 1 1
11 13 1 1
5 1 1 1
11 1 1 1
5 15 1 1
11 15 1 1
1 5 1 1
15 5 1 1
1 11 1 1
15 11 1 1
2 2 3 1
14 2 3 1
2 14 3 1
14 14 3 1
1 1 5 1
15 1 5 1
1 15 5 1
15 15 5 1
1 1 7 1
15 1 7 1
1 15 7 1
15 15 7 1
1 1 13 1
15 1 13 1
1 15 13 1
15 15 13 1
2 6 17 2
14 6 17 2
2 10 17 2
14 10 17 2
2 2 17 2
14 2 17 2
2 14 17 2
14 14 17 2
6 2 17 2
10 2 17 2
6 14 17 2
10 14 17 2
5 5 18 1
11 5 18 1
5 11 18 1
11 11 18 1
5 5 16 1
11 5 16 1
5 11 16 1
11 11 16 1
3 3 14 1
13 3 14 1
3 13 14 1
13 13 14 1
1 3 13 1
15 3 13 1
1 13 13 1
15 13 13 1
3 1 13 1
13 1 13 1
3 15 13 1
13 15 13 1
1 1 11 1
15 1 11 1
1 15 11 1
15 15 11 1
1 3 11 1
15 3 11 1
1 13 11 1
15 13 11 1
3 1 11 1
13 1 11 1
3 15 11 1
13 15 11 1
5 1 13 1
11 1 13 1
5 15 13 1
11 15 13 1
1 5 13 1
15 5 13 1
1 11 13 1
15 11 13 1
1 1 9 1
15 1 9 1
1 15 9 1
15 15 9 1
4 1 15 1
12 1 15 1
4 15 15 1
12 15 15 1
1 4 15 1
15 4 15 1
1 12 15 1
15 12 15 1
7 1 14 1
9 1 14 1
7 15 14 1
9 15 14 1
1 7 14 1
15 7 14 1
1 9 14 1
15 9 14 1
2 2 16 1
14 2 16 1
2 14 16 1
14 14 16 1
2 2 18 1
14 2 18 1
2 14 18 1
14 14 18 1
2 6 16 1
14 6 16 1
2 10 16 1
14 10 16 1
2 6 18 1
14 6 18 1
2 10 18 1
14 10 18 1
6 2 16 1
10 2 16 1
6 14 16 1
10 14 16 1
6 2 18 1
10 2 18 1
6 14 18 1
10 14 18 1
5 7 18 1
11 7 18 1
5 9 18 1
11 9 18 1
7 5 18 1
9 5 18 1
7 11 18 1
9 11 18 1
7 7 18 1
9 7 18 1
7 9 18 1
9 9 18 1];
disp(length(Points)+" spheres loaded");

% iterate over 187*186 combinations and check if any spheres overlap in 3D
iOverlaps=0;
iTotalVol=0;
for i = 1:length(Points)
    iTotalVol = iTotalVol + Points(i,4)^3;
    for j = 1:length(Points)
        if (i>j) && (sphereintersectvolume(norm(Points(i,1:3)-Points(j,1:3)),Points(i,4),Points(j,4)) ~= 0)
            iOverlaps = iOverlaps + 1;
            disp(["Sphere "+i+":"+Points(i,1)+"-"+Points(i,2)+"-"+Points(i,3)+"-"+Points(i,4)+", Sphere "+j+":"+Points(j,1)+"-"+Points(j,2)+"-"+Points(j,3)+"-"+Points(j,4)+", Distance: "+norm(Points(i,1:3)-Points(j,1:3))+", Intersecting Volume: "+sphereintersectvolume(norm(Points(i,1:3)-Points(j,1:3)),Points(i,4),Points(j,4))]);
        end
    end
end
disp(iOverlaps+" spheres found to intersect in 3D");

%find out the shortest distance of each potential point in the grid from each of the 187 Spheres
maxGap = 0;iPoints=0;
for i = 1:15
    for j = 1:15
        for k = 1:18
            sDist = 8;
            used = false;
            iPoints = iPoints +1;
            for m = 1:length(Points)
                if (i==Points(m,1) && j==Points(m,2) && k==Points(m,3))
                    used = true;
                elseif (abs(norm([i j k]-Points(m,1:3))-Points(m,4)) < sDist) 
                    sDist = abs(norm([i j k]-Points(m,1:3))-Points(m,4));
                end
            end
% check if the shortest distance is >1 unit
            if (not(used)) && (sDist >= 0)
                maxGap = max(maxGap,sDist);
%                disp(["Point:"+i+"~"+j+"~"+k+", Distance: "+sDist]);
            end
        end
    end
end
disp(iPoints+" points scanned");
disp("maximum internal/external gap is " + maxGap);
disp("total V is "+iTotalVol);

function V = sphereintersectvolume(D,R1,R2)
% V = sphereintersectvolume(D,R1,R2)
% arguments: (input) 
% D - distance between centers of the two spheres
% R1, R2 - radii of the pair of spheres 
if any([D,R1,R2] < 0)
  error('all radii and distances should in general be non-negative numbers.')
end
% we need to know which radius is smaller of the two
Rmin = min(R1,R2);
Rmax = max(R1,R2);
if D > (Rmin + Rmax)
  V = 0;
elseif Rmax >= (D + Rmin)
  % this triggers if one sphere lies entirely inside the other.
  % then the volume returned is just the volume of the smaller sphere.
  % let's change it to return 0 instead!
  V = 0;
else
  % There is a non-trivial intersection
  X = (D^2 - Rmax^2 + Rmin^2)/(2*D);
  V = spheresegmentvolume([X,Rmin],3,Rmin) + spheresegmentvolume([D - X,Rmax],3,Rmax);
end
end

function V = spheresegmentvolume(t,n,radius)
% spheresegmentvolume: n-d volume of a sphere cap or within any band defined by parallel planes
% usage: V = spheresegmentvolume(t,n,radius)
%
% Integrates a (hyper)sphere cap, or any segment of a sphere defined by
% a pair of parallel planar slices through the sphere.
% 
% arguments: (input)
%  t - sphere cap limit(s). t must be a vector of length 2
%      Each element of t must normally lie in the closed
%      interval  [-radius,radius].
% 
%      t is a vector that defines the spacing of the parallel
%      slicing planes along any of the axes. Of course, since
%      a sphere is fully rotationally symmetric, the actual
%      planes can lie along any axis. This allows us to define
%      the planes by a simple vector of length 2, denoting their
%      relative displacement from each other.
%
%      A sphere cap would be defined by the vector [t1,radius].
%      
%      If t is empty, or not provided, then the complete sphere
%      volume is computed.
%
%      if t(1) > t(2), then the computed volume will be negative.
%
%      Default value: t = [-radius,radius]. This will compute
%      the complete enclosed volume in the sphere.
% 
%  n - (OPTIONAL) - (scalar, positive numeric integer) dimension
%      of the n-sphere.
%
%      n must be an integer, greater than zero.
%
%      Note: When n is very large (on the order of 344) an underflow
%      will result, causing the computed volume to be returned
%      as zero.
%
%      Default value: n = 3
%
%  radius - (OPTIONAL) - (scalar numeric) radius of the hyper-sphere
%
%      radius must be non-negative
% 
%      Default value: radius = 1
%
%
% Arguments: (output)
%  V - n-dimensional volume of the indicated sphere cap or band.
%
%      The computed volume is exact, to within the floating
%      point accuracy available in double precision.
%
%
% Example:
% %  Compute the volume of a unit, complete sphere
% %  In 2-d, the "volume" is pi.
%
%  V = spheresegmentvolume([],2)
% % V =
% %       3.1416
%
%  V - pi
% % ans =
% %       0
%
% % In 3-d, compute the volume of an exact unit hemisphere,
% % the volume is pi*2/3, approximately 2.0944
%
%  V = spheresegmentvolume([0,1],3)
% % V =
% %       2.0944
%
%  V - 2*pi/3
% % ans =
% %      -4.4409e-16
%   
% % In 4-d, compute the volume of a band around a hyper-sphere
% % of radius 2, with the band running from -1 to +1.
%
%  V = spheresegmentvolume([-1,1],4,2)
% % V =
% %      58.967
%
% % In 50 dimensions, compute the volume inside a unit hemi-sphere cap
%  
%  V = spheresegmentvolume([0,1],50)
% % V =
% %       8.6511e-14
%
% See also: 
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 2.0
% Release date: 2/05/08
% apply defaults as indicated
if (nargin<3) || isempty(radius)
  radius = 1;
elseif length(radius)>1
  error('radius must be a scalar if provided')
elseif radius < 0
  error('radius must be non-negative')
elseif radius == 0
  % a zero radius implies a zero volume always.
  V = 0;
  return
end
if (nargin<2) || isempty(n)
  n = 3;
elseif (length(n)>1) || (n<1) || (n~=floor(n))
  error('n must be a positive scalar integer if provided')
end
% if t is empty, then the complete sphere volume will be computed
if (nargin<1) || isempty(t)
  t = [];
end
t = t(:);
if (length(t)>2) || (length(t) == 1)
  error('t must be a vector of length 2 if provided')
end
% we are done if t is empty
if isempty(t)
  % a complete sphere
  V = (radius^n)*2*(pi^(n/2))/n/gamma(n/2);
  
elseif n==1
  V = diff(t);
else
  % its only a segment from a sphere
  
  % compute the complete sphere volume in one lower dimension
  nm1 = n-1;
%  V = 2*(pi^(nm1/2))/nm1/gamma(nm1/2);
  V = (pi^(nm1/2))/gamma(nm1/2 + 1);
  
  % non-dimensionalize the problem
  t = t/radius;
  
  % ensure that both endpoints are in the nondimensional
  % interval [-1,1]
  t = min(1,max(-1,t));
  
  % this problem reduces to an integral of sin(acos(t))^(n-1), done by
  % integration by parts. Transformation of variables, by s = acos(t)
  % turns that into an integral of -sin(s)^n
  
  % transform the limits
  s = acos(t);
  V = V*(radius^n)*-sineintegral(n,s);
end  
end

function si = sineintegral(n,s)
% integral of sin(s)^n, with limits of integration [s(1),s(2)].
% done recursively.
if n == 1
  % special case when n == 1 to end the recurrence.
  si = diff(-cos(s));
elseif n==2
  % special case when n == 2
  si = diff(s/2 - sin(2*s)/4);
else
  % n must be greater than 2
  si = diff(-cos(s).*sin(s).^(n-1)/n) + ((n-1)/n)*sineintegral(n-2,s);
end
end