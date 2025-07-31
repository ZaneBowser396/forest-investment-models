function carbon_data = carbon_reader(filename, bbox)
fid = fopen(filename, 'r');
header_info = textscan(fid, '%s %f', 6);

ncols = header_info{2}(1);
nrows = header_info{2}(2);
xllcorner = header_info{2}(3);
yllcorner = header_info{2}(4);
cellsize = header_info{2}(5);
nodata_value = header_info{2}(6);

data = fscanf(fid, '%f', [ncols, nrows]);
data = data';
fclose(fid);
data(data == nodata_value) = NaN;

lon_asc = xllcorner + (0:ncols-1) * cellsize;
lat_asc = yllcorner + (0:nrows-1) * cellsize;
lat_asc = flip(lat_asc);

lon_indices = (lon_asc >= bbox(2,1)) & (lon_asc <= bbox(1,1));
lat_indices = (lat_asc >= bbox(2,2)) & (lat_asc <= bbox(1,2));

filtered_data = data(lat_indices, lon_indices);
filtered_lon = lon_asc(lon_indices);
filtered_lat = lat_asc(lat_indices);

carbon_data = struct('Data', filtered_data, 'Lon', filtered_lon, 'Lat', filtered_lat, 'Cellsize', cellsize);

end

