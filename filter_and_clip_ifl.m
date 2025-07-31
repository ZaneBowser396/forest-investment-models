%% Function to filter and clip IFL polygons based on bounding box and intersection with PNG border
function filtered_data = filter_and_clip_ifl(ifl_data, png_data, bound_box)
    filtered_data = [];
    png_poly = polyshape(png_data(1).X, png_data(1).Y);
    
    for i = 1:length(ifl_data)
        if isfield(ifl_data(i), 'X') && isfield(ifl_data(i), 'Y')
            X = ifl_data(i).X;
            Y = ifl_data(i).Y;
        
            if any(X >= bound_box(2,1) & X <= bound_box(1,1)) && ...
               any(Y >= bound_box(2,2) & Y <= bound_box(1,2))
                
                ifl_poly = polyshape(X, Y);
                clipped_poly = intersect(ifl_poly, png_poly);
                
                if ~isempty(clipped_poly.Vertices)
                    filtered_data = [filtered_data; struct('poly', clipped_poly, 'IFL_ID', ifl_data(i).IFL_ID, 'AREA2020', ifl_data(i).AREA2020)];
                end
            end
        end
    end
end