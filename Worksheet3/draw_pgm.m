%
% drawpicture.m - Create a PGM image file using pgmwrite.m
%
 

% % KARMAN: hardcoded approach
% % Create a matrix of the desired image dimensions
% % In this case, height is xx and width is xx
% % Initialize p values to a +ve number for fluid region
% height = 20;
% width = 100;
% p = 255*ones(height,width);
 
% p(3*height/5-1:3*height/5, 2*height/5+1) = 0;
% p(3*height/5-2:3*height/5, 2*height/5+2) = 0;
% p(2*height/5+1:2*height/5+3, 2*height/5+3) = 0;
% p(2*height/5+1:2*height/5+2, 2*height/5+4) = 0;
% % End of KARMAN


% % SHEAR
% % Create a matrix of the desired image dimensions
% % In this case, height is xx and width is xx
% % Initialize p values to a +ve number for fluid region
% height = 20;
% width = 100;
% p = 255*ones(height,width);
 
% % End of SHEAR


% % STEP
% % Create a matrix of the desired image dimensions
% % In this case, height is xx and width is xx
% % Initialize p values to a +ve number for fluid region
% height = 20;
% width = 100;
% p = 255*ones(height,width);
 
% % Calculate the pixel values using a nested while loop
% y = 1;
% while y <= height% % Create a matrix of the desired image dimensions
% % In this case, height is xx and width is xx
% % Initialize p values to a +ve number for fluid region

%     x = 1;
%     while x <= width
%         if y > height/2 && x <= height/2
%             p(y,x) = 0;
%         end
%         x = x + 1;
%     end
%     y = y + 1;
% end
% % End of STEP


% % CAVITY
% % Create a matrix of the desired image dimensions
% % In this case, height is xx and width is xx
% % Initialize p values to a +ve number for fluid region
height = 40; %20;
width = 50; %100;
p = 255*ones(height,width);
 
% % End of CAVITY

% % Write pixels to a PGM image file
pgmwrite(p);