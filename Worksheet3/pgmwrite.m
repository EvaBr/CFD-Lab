% pgmwrite.m - MATLAB function to write a matrix
% of pixel values to a greyscale PGM image file.
%
% Written by Ted Burke - last updated 12-3-2015
%
 
function pgmwrite(pixels)
    % Check matrix dimensions which will determine
    % the width and height of the image
    s = size(pixels);
    height = s(1);
    width = s(2);
     
    % Open a file for writing
    fid = fopen('image.pgm','w');
     
    % Write the PGM image header to the file
    fprintf(fid,'P2\n');
    fprintf(fid,'%d %d\n', width, height);
    fprintf(fid,'255\n');
     
    % Write the pixel values from the matrix into
    % the file
    y = 1;
    while y <= height
        x = 1;
        while x <= width
            fprintf(fid, '%03d ', pixels(y,x));
            x = x + 1;
        end
        fprintf(fid, '\n');
        y = y + 1;
    end
     
    % Close the file
    fclose(fid)
end