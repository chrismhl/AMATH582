%vidFrames is the .dat file corresponding to the recording
%guess is the initial guess of the flashlight in [x y] from frame 1
%width is the search width of the tracking algorithm
%returns vectors xp and yp which are the positions of the paintcan with
%   each frame
function [xp,yp] = track_paintcan(vidFrames,guess,width)
dim = size(vidFrames);
L = dim(4);     %length of video in frames
xp = zeros(L,1);
yp = zeros(L,1);
maxloc = [0,0];
maxval = 0;

init_guess = guess;
xp(1) = init_guess(1);
yp(1) = init_guess(2);

%center the first search window to the initial point
center = init_guess;

for f = 2:L
    frame = double(vidFrames(:,:,:,f));
     for x = center(1)-width(1):center(1)+width(2)
        for y = center(2)-width(3):center(2)+width(4)
            point = frame(y,x,:);
            colorsum = point(1,1,1)+point(1,1,2)+point(1,1,3);
            if colorsum > maxval
                maxloc = [x,y];
                maxval = colorsum;
            end   
        end
    end
    xp(f) = maxloc(1);
    yp(f) = maxloc(2);
    maxloc = [0,0];
    maxval = 0;
    
    center(1) = xp(f);
    center(2) = yp(f);
end
