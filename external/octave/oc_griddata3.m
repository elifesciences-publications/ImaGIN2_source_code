%# -*- texinfo -*-
%# @deftypefn {Function File} {@var{vi} =} griddata3 (@var{x}, @var{y}, @var{z}, @var{v} @var{xi}, @var{yi}, @var{zi}, @var{method}, @var{options})
%# 
%# Generate a regular mesh from irregular data using interpolation.
%# The function is defined by @code{@var{y} = f (@var{x},@var{y},@var{z})}.
%# The interpolation points are all @var{xi}.  
%#
%# The interpolation method can be @code{'nearest'} or @code{'linear'}.
%# If method is omitted it defaults to @code{'linear'}.
%# @seealso{griddata, delaunayn}
%# @end deftypefn

%# Author: David Bateman <dbateman@free.fr>

function vi = oc_griddata3 (x, y, z, v, xi, yi, zi, method, varargin)
        
if (nargin < 7)
    print_usage ;
end

if (~all (size (x) == size (y) & size (x) == size(z) & size(x) == size (v)))
    error ('griddata3: x, y, z, and v must be vectors of same length');
end

%# meshgrid xi, yi and zi if they are vectors unless they
%# are vectors of the same length
if (isvector (xi) && isvector (yi) && isvector (zi) && (numel (xi) ~= numel (yi) || numel (xi) ~= numel (zi)))
    [xi, yi, zi] = meshgrid(xi, yi, zi);
end

if (any (size(xi) ~= size(yi)) || any (size(xi) ~= size(zi)))
    error ('griddata3: xi, yi and zi must be vectors or matrices of same size');
end

vi = oc_griddatan([x(:), y(:), z(:)], v(:), [xi(:), yi(:), zi(:)], varargin{:});
vi = reshape(vi, size (xi));


  