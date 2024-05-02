function [domain, lat, lon] = domain_subset(lat, lon, ...
            min_lat, max_lat, min_lon, max_lon)

% Get a subset of grid cells only
% Gabrielle De Lannoy - 13 March 2023
%--------------------------------------------------------------------------

if (~isnan(min_lon) && ~isnan(min_lat) && ...
    ~isnan(min_lat) && ~isnan(max_lat))
    domain  = find(lon>=min_lon & lon<=max_lon & ...
               lat>=min_lat & lat<=max_lat);
    lat     = lat(domain);
    lon     = lon(domain);
else
    domain = [1:length(lat)];
end

end
%==============EOF================================
