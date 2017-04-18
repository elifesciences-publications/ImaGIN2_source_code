%# -*- texinfo -*-
%# @deftypefn {Function File} {@var{yi} =} griddatan (@var{x}, @var{y}, @var{xi}, @var{method}, @var{options})
%# 
%# Generate a regular mesh from irregular data using interpolation.
%# The function is defined by @code{@var{y} = f (@var{x})}.
%# The interpolation points are all @var{xi}.  
%#
%# The interpolation method can be @code{'nearest'} or @code{'linear'}.
%# If method is omitted it defaults to @code{'linear'}.
%# @seealso{griddata, delaunayn}
%# @end deftypefn

%# Author: David Bateman <dbateman@free.fr>

function yi = oc_griddatan (x, y, xi, method, varargin)

if (nargin == 3)
    method = 'linear';
end
if (nargin < 3)
    print_usage ;
end

if (ischar (method))
    method = lower (method);
end

[m, n] = size (x);
[mi, ni] = size (xi);

if (n ~= ni || size (y, 1) ~= m || size (y, 2) ~= 1)
    error ('griddatan: dimensional mismatch');
end

%# triangulate data
%# tri = delaunayn(x, varargin{:});
tri = delaunayn(x);

yi = NaN(mi, 1);

if (strcmp (method, 'nearest'))
    %# search index of nearest point
    idx = dsearchn (x, tri, xi);
    valid = ~isnan (idx);
    yi(valid) = y(idx(valid));
    
elseif (strcmp (method, 'linear'))
    %# search for every point the enclosing triangle
    [tri_list, bary_list] = tsearchn (x, tri, xi);
    
    %# only keep the points within triangles.
    valid = ~isnan (tri_list);
    tri_list = tri_list(~isnan (tri_list));
    bary_list = bary_list(~isnan (tri_list), :);
    nr_t = size(tri_list,1);
    
    %# Use barycentric coordinate of point to calculate yi
    yi(valid) = sum (y(tri(tri_list,:)) .* bary_list, 2);
    
else
    error ('griddatan: unknown interpolation method');
end



