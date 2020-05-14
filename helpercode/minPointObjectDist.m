function [x,y,min_dist] = minPointObjectDist(row1,col1,comp2)
    [rows2,cols2] = ind2sub(size(comp2),find(bwperim(comp2)));

    [min_dist,min_ind2] = min(sqrt((rows2 - row1).^2 ...
        + (cols2 - col1).^2));
    
    [x,y] = bresenham(rows2(min_ind2), cols2(min_ind2), ...
        row1, col1);
end