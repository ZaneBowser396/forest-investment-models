%% Setup
clear, clc
warning('off', 'MATLAB:polyshape:repairedBySimplify');
warning('off', 'MATLAB:polyshape:boolOperationFailed');

%% Read in datasets

%% PNG Boundary
disp('Setting initial PNG Data');
bound_box = [159.4929166 -0.872038700000000; 140.841655300000 -11.6547069000000];
PNG_boundaries = shaperead('geoBoundaries-PNG-ADM0-all/geoBoundaries-PNG-ADM0.shp');
% https://media.githubusercontent.com/media/wmgeolab/geoBoundaries/9469f09592ced973a3448cf66b6100b741b64c0d/releaseData/gbOpen/PNG/ADM0/geoBoundaries-PNG-ADM0-all.zip
PNG_border = polyshape(PNG_boundaries.X, PNG_boundaries.Y);

%% IFL Datasets
disp('Reading IFL data');
ifls_2020 = shaperead('IFL_2020/ifl_2020.shp');
ifls_2016 = shaperead('IFL_2016/ifl_2016.shp');
ifls_2013 = shaperead('IFL_2013/ifl_2013.shp');

disp('Filtering IFLs to only PNG');
filtered_2020 = filter_and_clip_ifl(ifls_2020, PNG_boundaries, bound_box);
filtered_2016 = filter_and_clip_ifl(ifls_2016, PNG_boundaries, bound_box);
filtered_2013 = filter_and_clip_ifl(ifls_2013, PNG_boundaries, bound_box);

disp('Sorting IFLs by ID')
sorted_2020 = sort_ifl(filtered_2020);
sorted_2013 = sort_ifl(filtered_2013);
sorted_2016 = sort_ifl(filtered_2016);

%% Safety Datasets

disp('Reading Safety Data');
lat = -12:0.025:1;
lon = 140:0.025:157;
[X,Y] = meshgrid(lon,lat);

Safety = readmatrix('PNG_Safety.csv',Range=7);

%% Accessibility Datasets

disp('Reading Accessibility Data');
roads = shaperead("groads-v1-oceania-east-shp/gROADS-v1-oceania-east.shp", "UseGeoCoords", true);
airstrip = shaperead("MAF_PNG_Rural_Airstrip_Status_20222/MAF_PNG_Rural_Airstrip_Status_20222.shp", "UseGeoCoords", true);
port_moresby = [147.16641, -9.42257];

% Select only the roads and airstrips inside PNG;
filtered_roads = roads(arrayfun(@(r) ...
    any(r.Lon >= bound_box(2,1) & r.Lon <= bound_box(1,1) & ...
        r.Lat >= bound_box(2,2) & r.Lat <= bound_box(1,2) & ...
        (~(r.Lon < 145) | (r.Lat > -10))), roads));

filtered_airstrips = airstrip(arrayfun(@(a) ...
    any(a.Lon >= bound_box(2,1) & a.Lon <= bound_box(1,1) & ...
        a.Lat >= bound_box(2,2) & a.Lat <= bound_box(1,2) & ...
        a.StatusFor == "Open for MAF Operations" & ...
        (~(a.Lon < 145) | (a.Lat > -10))), airstrip));

%% Carbon Dataset

disp('Reading Carbon Data')
above_carbon_data = carbon_reader('HDCAbove01.asc', bound_box);

%% Species Dataset

disp('Reading Species Data');
species_data = readtable('SpeciesRanges.csv');


%% Calculate IFL Data

%% Calculate Area

disp('Storing the area of each IFL');
ids = cell(length(sorted_2020), 1);
A_2013 = NaN(length(sorted_2020),1);
A_2016 = NaN(length(sorted_2020),1);
A_2020 = zeros(length(sorted_2020),1);

for i = 1:length(sorted_2020)
    % Get the ID of the current IFL
    ids{i} = sorted_2020(i).IFL_ID;
    A_2020(i) = sorted_2020(i).AREA2020;

    % Find the corresponding Area for 2013 and 2016 IFLS
    % Problems with 125_1 vs 125, to be fixed
    idx_2013 = find(strcmp({sorted_2013.IFL_ID}, ids{i}), 1);
    if ~isempty(idx_2013)
        A_2013(i) = sorted_2013(idx_2013).AREA2020;
    end

    idx_2016 = find(strcmp({sorted_2016.IFL_ID}, ids{i}), 1);
    if ~isempty(idx_2016)
        A_2016(i) = sorted_2016(idx_2016).AREA2020;
    end
end

%% Calculate Safety

disp('Calculating average safety');
avg_safety = zeros(length(sorted_2020),1);
Safety_no_nan = Safety;
Safety_no_nan(isnan(Safety)) = 0;

Safety_lat_range = linspace(lat(1), lat(end), size(Safety,1));
Safety_lon_range = linspace(lon(1), lon(end), size(Safety,2));

for i=1:length(sorted_2020)
    ifl_poly = sorted_2020(i).poly;

    % Get the safety coords that are in the IFL
    in = inpolygon(X, Y, ifl_poly.Vertices(:,1), ifl_poly.Vertices(:,2));

    [lat_idx, lon_idx] = find(in);

    % Clip the lat/long coords to the nearest safety value
    nearest_lat_idx = interp1(Safety_lat_range, 1:size(Safety,1), Y(in), 'nearest', 'extrap');
    nearest_lon_idx = interp1(Safety_lon_range, 1:size(Safety,2), X(in), 'nearest', 'extrap');
    flipped_lat_idx = size(Safety,1) - nearest_lat_idx + 1;

    % Get the safety values that match coords in the IFL
    safety_values = Safety_no_nan(sub2ind(size(Safety_no_nan), flipped_lat_idx, nearest_lon_idx));

    total_safety = sum(safety_values);
    count = numel(safety_values);
    
    % Average the safety values, NaN if there are none (should not happen
    % unless scale is too big)
    if count > 0
        avg_safety(i) = total_safety / count;
    else
        avg_safety(i) = NaN;
    end
end

%% Calculate Carbon Data

disp('Calculating average carbon');
above_carbon = zeros(length(sorted_2020), 1);
[carbonX, carbonY] = meshgrid(above_carbon_data.Lon, above_carbon_data.Lat);
data_scale = 1/100; % Proper scaling of units

for i=1:length(sorted_2020)
    poly = sorted_2020(i).poly;

    % Get the carbon indexes inside the IFL
    inside = inpolygon(carbonX, carbonY, poly.Vertices(:,1), poly.Vertices(:,2));

    % Pick out the carbon values inside the IFL
    [lat_idx, lon_idx] = find(inside);
    carbon_values_inside = above_carbon_data.Data(sub2ind(size(above_carbon_data.Data), lat_idx, lon_idx));

    total_carbon = sum(carbon_values_inside);
    count = numel(carbon_values_inside);

    % Average scaled carbon
    if count > 0
        above_carbon(i) = total_carbon / count * data_scale;
    else
        above_carbon(i) = NaN;
    end
end

%% Calculate Travel Time

disp('Calculating travel time');
normalTravel = zeros(length(sorted_2020), 1);
airTravel = zeros(length(sorted_2020), 1);
travelTime = cell(length(sorted_2020), 2);

for i=1:length(sorted_2020)
    ifl_poly = sorted_2020(i).poly;

    % Find the nearest point on the IFL to port moresby, and the distance
    [vertexID] = nearestvertex(ifl_poly, port_moresby);
    distance = haversine(port_moresby(2), port_moresby(1), ifl_poly.Vertices(vertexID,2), ifl_poly.Vertices(vertexID,1));

    % Split into evenly spaced steps
    num_steps = 300; % best accuracy to time i could be bothered to find
    lats_path = linspace(port_moresby(2), ifl_poly.Vertices(vertexID,2), num_steps);
    lons_path = linspace(port_moresby(1), ifl_poly.Vertices(vertexID,1), num_steps);

    total_travel_time = 0;

    % Approximate each step with either inside PNG or outside, ignore roads
    % for now
    for j = 1:num_steps-1
        segment_distance = haversine(lats_path(j), lons_path(j), lats_path(j+1), lons_path(j+1));
        inside_png = isinterior(PNG_border, lons_path(j), lats_path(j));
        
        if inside_png
            speed = 4;
        else
            speed = 20;
        end
        
        travel_time = segment_distance / speed;
        total_travel_time = total_travel_time + travel_time;
    end
    
    normalTravel(i) = total_travel_time;


    % Calculate comparison time by flying to nearest airstrip then driving
    % to IFL
    min_airstrip_distance = Inf;
    nearest_airstrip = [];

    for k = 1:length(filtered_airstrips)
        airstrip_lat = filtered_airstrips(k).Lat;
        airstrip_lon = filtered_airstrips(k).Lon;

        [vertexID] = nearestvertex(ifl_poly, [airstrip_lon, airstrip_lat]);
        airstrip_distance = haversine(airstrip_lat, airstrip_lon, ifl_poly.Vertices(vertexID, 2), ifl_poly.Vertices(vertexID,1));

        if airstrip_distance < min_airstrip_distance
            min_airstrip_distance = airstrip_distance;
            nearest_airstrip = [airstrip_lat, airstrip_lon, vertexID];
        end
    end

    % Two days to fly to airstrip
    total_flight_time = 48;
    lats_path = linspace(nearest_airstrip(1), ifl_poly.Vertices(nearest_airstrip(3), 2), num_steps);
    lons_path = linspace(nearest_airstrip(2), ifl_poly.Vertices(nearest_airstrip(3), 1), num_steps);

    % Step to IFL from airstrip and approximate time
    for j=1:num_steps-1
        segment_distance = haversine(lats_path(j), lons_path(j), lats_path(j+1), lons_path(j+1));
        inside_png = isinterior(PNG_border, lons_path(j), lats_path(j));

        if inside_png
            speed = 4;
        else
            speed = 20;
        end
        
        flight_time = segment_distance / speed;
        total_flight_time = total_flight_time + flight_time;
    end

    airTravel(i) = total_flight_time;

    % Use the fastest time
    if normalTravel(i) >= airTravel(i)
        travelTime{i,1} = airTravel(i);
        travelTime{i,2} = "Flight";
    else
        travelTime{i,1} = normalTravel(i);
        travelTime{i,2} = "No Flight";
    end

end

%% Calculate Species Richness

disp('Calculating species richness')

species_count = zeros(length(sorted_2020), 1);

columnNames = species_data.Properties.VariableNames;

% Seperate each species out to seperate names
for i=3:numel(columnNames)
    columnData = species_data.(columnNames{i});
    columnType = class(columnData);
    if strcmp(columnType, 'cell')
        numericData = cellfun(@str2double, columnData);
        species_data.(columnNames{i}) = numericData;
    end
end

speciesPolyshapes = containers.Map();
speciesNames = {};
speciesTypes = {};

% Get a single polyshape for each species, combining multiple ones
for i = 1:2:height(species_data)
    idx = (i+1)/2;

    currentSpeciesName = species_data.Var1{i};
    currentSpeciesType = species_data.Var2{i};

    longitudes = species_data{i, 3:end};
    latitudes = species_data{i+1, 3:end};

    % Get lat/long coords up until first NaN or end
    nanIdx = find(isnan(longitudes), 1);
    if isempty(nanIdx)
        nanIdx = length(longitudes) + 1;
    end

    % Get outer boundary of polyshape up until first NaN
    outerBoundary = polyshape(longitudes(1:nanIdx-1), latitudes(1:nanIdx-1));

    % Iterate over remaining code, adding inner boundaries up to every NaN
    % until no remaining coords to be used
    innerBoundaries = {};
    if nanIdx <= length(longitudes)
        remainingLongitudes = longitudes(nanIdx+1:end);
        remainingLatitudes = latitudes(nanIdx+1:end);
        
        while ~isempty(remainingLongitudes) && ~all(isnan(remainingLongitudes))
            nanIdx = find(isnan(remainingLongitudes), 1);
            if isempty(nanIdx)
                nanIdx = length(remainingLongitudes) + 1;
            elseif all(isnan(remainingLongitudes(1:nanIdx-1)))
                break;
            end
            
            innerBoundary = polyshape(remainingLongitudes(1:nanIdx-1), remainingLatitudes(1:nanIdx-1));
            innerBoundaries{end+1} = innerBoundary;
            
            if nanIdx < length(remainingLongitudes)
                remainingLongitudes = remainingLongitudes(nanIdx+1:end);
                remainingLatitudes = remainingLatitudes(nanIdx+1:end);
            else
                remainingLongitudes = [];
                remainingLatitudes = [];
            end
        end
    end

    % Add every inner boundary to the polyshape
    for j=1:length(innerBoundaries)
        outerBoundary = addboundary(outerBoundary, innerBoundaries{j}.Vertices(:,1), innerBoundaries{j}.Vertices(:,2));
    end

    % If any of the species poly is inside PNG, add it to the map,
    % combining if exists or making a new entry if it doesnt
    in = intersect(outerBoundary, PNG_border);
    if area(in) > 0
        if isKey(speciesPolyshapes, currentSpeciesName)
            speciesPolyshapes(currentSpeciesName) = union(speciesPolyshapes(currentSpeciesName), outerBoundary);
        else
            speciesPolyshapes(currentSpeciesName) = outerBoundary;
            speciesNames{end+1} = currentSpeciesName;
            speciesTypes{end+1} = currentSpeciesType;
        end
    end
end

% For every IFL, for every species, if the species is in the poly add one
% to the count
for i=1:length(sorted_2020)
    poly = sorted_2020(i).poly;
    uniqueSpeciesInLandscape = containers.Map();

    speciesKeys = keys(speciesPolyshapes);
    for j=1:length(speciesKeys)
        currentSpecies = speciesKeys{j};
        speciesShape = speciesPolyshapes(currentSpecies);

        in = intersect(poly, speciesShape);
        if ~isempty(in.Vertices)
            if ~isKey(uniqueSpeciesInLandscape, currentSpecies)
                uniqueSpeciesInLandscape(currentSpecies) = true;
            end
        end
    end

    species_count(i) = length(uniqueSpeciesInLandscape);
end

%% Static Ranking

%% Parameter Calculation

% Calculate loss rate - assume difference over four years is constant rate
delta = (A_2016 - A_2020) ./ (A_2016 * 4);

% Calculate probability of local support
p = 1 - exp(-25 * delta .* above_carbon);

% Calculate binomial probability of success based on less than two 
% incidents within 20 visits, extrapolated from safety score
q = zeros(length(avg_safety),1);
n_visits = 20;

for i=1:length(avg_safety)
    qi = interp1([0 1 2 3 4 5 6], [0 0.01 0.02 0.04 0.1 0.33 0.5], avg_safety(i));
    
    binomial_coeff = @(n,k) factorial(n)/(factorial(k)*factorial(n-k));

    Pr_0 = binomial_coeff(n_visits, 0) * (qi^0) * (1-qi)^(n_visits-0);
    Pr_1 = binomial_coeff(n_visits, 1) * (qi^1) * (1-qi)^(n_visits-1);
    
    q(i) = Pr_0 + Pr_1;
end

%% Static Score Calculation
E = p .* q .* species_count / travel_time;

%% Plot Data

% Accessability Plot
figure
hold on;
plot(PNG_border)
xlim([bound_box(2,1), bound_box(1,1)])
ylim([bound_box(2,2), bound_box(1,2)])
for i=1:length(filtered_roads)
    plot(filtered_roads(i).Lon, filtered_roads(i).Lat, 'r');
end
for i = 1:length(filtered_airstrips)
    plot(filtered_airstrips(i).Lon, filtered_airstrips(i).Lat, 'o', 'MarkerSize', 2, 'MarkerEdgeColor', 'magenta')
end
xlabel('Longitude');
ylabel('Latitude');
title('Roads and Airstrips of PNG');
% legend({'Land', 'Roads', 'Airstrips'})
plot(port_moresby(1), port_moresby(2), 'o', 'MarkerSize', 15, 'MarkerEdgeColor', 'green');
hold off;

% IFLs
figure
hold on;
plot(PNG_border);

for i = 1:length(sorted_2020)
    poly = sorted_2020(i).poly;
    plot(poly, 'FaceColor', 'green');
end
xlim([bound_box(2,1), bound_box(1,1)])
ylim([bound_box(2,2), bound_box(1,2)])
title('Intact Forest Landscapes');
xlabel('Longitude');
ylabel('Latitude');
hold off;

% Safety Heatmap
figure
hold on;
surf(lon(1:2:end),lat(1:2:end),flipud(Safety));
colormap('spring');
colorbar
shading interp
view(2)
xlabel('Longitude');
ylabel('Latitude')
title('Expert Safety Scores for PNG');
hold off;

% Species Areas
figure
hold on;
plot(PNG_border)
xlim([bound_box(2,1), bound_box(1,1)])
ylim([bound_box(2,2), bound_box(1,2)])
for i = 1:length(speciesPolyshapes)
    poly = speciesPolyshapes(speciesNames{i});
    plot(poly);
end
xlabel('Longitude');
ylabel('Latitude')
title('Species Locations PNG');
hold off;

% Carbon Plot
fid = fopen('HDCAbove01.asc', 'r');
% Read the header information
header_info = textscan(fid, '%s %f', 6);  % First 6 lines contain header info

% Extract the header information
ncols = header_info{2}(1);
nrows = header_info{2}(2);
xllcorner = header_info{2}(3);  % Longitude of lower-left corner
yllcorner = header_info{2}(4);  % Latitude of lower-left corner
cellsize = header_info{2}(5);   % Cell size (resolution)
nodata_value = header_info{2}(6);

% Read the grid data
data = fscanf(fid, '%f', [ncols, nrows]);
data = data';

% Close the file
fclose(fid);

% Replace no data values with NaN
data(data == nodata_value) = NaN;
% Define bounding box
x_min = 140;
x_max = 158;
y_min = -12;
y_max = 0;

% Calculate grids
lon_asc = xllcorner + (0:ncols-1) * cellsize;
lat_asc = yllcorner + (0:nrows-1) * cellsize;
lat_asc = flip(lat_asc);

% Filter data 
lon_indices = (lon_asc >= x_min) & (lon_asc <= x_max);
lat_indices = (lat_asc >= y_min) & (lat_asc <= y_max);

filtered_data = data(lat_indices, lon_indices);
filtered_lon = lon_asc(lon_indices);
filtered_lat = lat_asc(lat_indices);

% Visualize filtered data
figure;
imagesc(filtered_lon, filtered_lat, filtered_data);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('Longitude');
ylabel('Latitude');
title('Filtered Stored Carbon Data within Bounding Box');


%% Save Data

T = table(ids, A_2013, A_2016, A_2020, avg_safety, travelTime(:,1), travelTime(:,2), species_count, above_carbon, delta, p, q, E, ...
    'VariableNames', {'IFL_ID', 'A_2013', 'A_2016', 'A_2020', 'Average Safety', 'Travel_Time', 'Travel Type', 'Species_Richness', 'Average_Carbon', 'Loss Rate', 'Local Support', 'Chance of success', 'Overall ranking value'});
filename = 'tabulated_ifl_data.xlsx';
writetable(T,filename,'Sheet',1);

%% Useful Functions
function numericData = convertCellData(cellData)
    if isempty(cellData) || (iscell(cellData) && isempty(cellData{1}))
        numericData = NaN;  
    else
        if iscell(cellData)
            numericData = str2double(cellData{1}); 
        else
            numericData = cellData;  
        end
    end
end

function d = haversine(lat1, lon1, lat2, lon2)
    R = 6371;
    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    lat2 = deg2rad(lat2);
    lon2 = deg2rad(lon2);
    
    dlat = lat2 - lat1;
    dlon = lon2 - lon1;
    
    a = sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    
    d = R * c;
end

