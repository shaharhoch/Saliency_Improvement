function [ out_image ] = interpolateImage( in_image, image_boarders, max_rand, extrapolation_boarders )
    %Define Random number macros 
    randNumPos = @() randi([0, max_rand]);
    randNum = @() randi([(-1)*max_rand, max_rand]);
    
    out_image = in_image; 
    
    for row=extrapolation_boarders(1,1):extrapolation_boarders(1,2)
        for col=extrapolation_boarders(2,1):extrapolation_boarders(2,2)
            if(isInsideBaorders([row,col], image_boarders) == true)
                continue;
            end
            
            %Find the row index of the replacement pixel
            if((row >= image_boarders(1,1)) && (row <= image_boarders(1,2)))
                row_replace = row + randNum(); 
            elseif(row < image_boarders(1,1))
                row_replace = image_boarders(1,1) + randNumPos();
            else
                row_replace = image_boarders(1,2) - randNumPos();
            end
            % Make sure the row is in the picture bounds
            row_replace = max([row_replace, image_boarders(1,1)]);
            row_replace = min([row_replace, image_boarders(1,2)]);
            
            %Find the column index of the replacement pixel
            if((col >= image_boarders(2,1)) && (col <= image_boarders(2,2)))
                col_replace = col + randNum(); 
            elseif(col < image_boarders(2,1))
                col_replace = image_boarders(2,1) + randNumPos();
            else
                col_replace = image_boarders(2,2) - randNumPos();
            end
            % Make sure the col is in the picture bounds
            col_replace = max([col_replace, image_boarders(2,1)]);
            col_replace = min([col_replace, image_boarders(2,2)]);
            
            %Replace the pixel
            out_image(row, col, :) = in_image(row_replace, col_replace, :);
        end
    end


end

function [ is_inside ] = isInsideBaorders( point, image_boarders )
    is_inside = true;
    
    if((point(1) < image_boarders(1,1)) || (point(1) > image_boarders(1,2)))
        is_inside = false;
        return;
    end
    
    if((point(2) < image_boarders(2,1)) || (point(2) > image_boarders(2,2)))
        is_inside = false;
        return;
    end
end