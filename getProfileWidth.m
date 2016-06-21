function [ width, x1, x2 ] = getProfileWidth( x, y, varargin )

options = struct( 'smoothspan', 1, ...
	'baseline', 0, ...
	'type', 'rms' );
validProperties = fieldnames( options );

nArgs = length( varargin );
if mod( nArgs, 2 )
    error('getProfileWidth needs propertyName/propertyValue pairs')
end

for pair = reshape( varargin, 2, [] )
   property = lower( pair{1} );
   if any( strcmp( property, validProperties ) )
       options.(property) = pair{2};
   else
       error( '%s is not a recognized property name', property )
   end
end

if isrow(x); x = x'; end
if isrow(y); y = y'; end
y = smooth( y, options.smoothspan ) - options.baseline;

switch options.type
	case 'rms'
		w = y/sum(y);
		xMean = sum(x.*w);
		width = sqrt( sum( (x-xMean).^2.*w ) );
		x1 = xMean - width;
		x2 = xMean + width;
	case 'fwhm'
		yMax = max(y);
		idx1 = find( y > 0.5 * yMax, 1, 'first' );
		idx2 = find( y > 0.5 * yMax, 1, 'last' );
		x1 = interp1( y( idx1-1:idx1 ), x( idx1-1:idx1 ), 0.5 * yMax );
		x2 = interp1( y( idx2:idx2+1 ), x( idx2:idx2+1 ), 0.5 * yMax );
		width = abs( x2-x1 );
	otherwise
		error( 'invalid type includes rms and fwhm' );
end

end

