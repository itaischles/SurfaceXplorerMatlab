function newmap = bluewhitered(m, lims)

bottom = [0 0 0.2];
botmiddle = [0 0.1 1];
middle = [1 1 1];
topmiddle = [1 0 0];
top = [0.4 0 0];
% Find middle
% % % lims = get(gca, 'CLim');
% Find ratio of negative to positive
if (lims(1) < 0) & (lims(2) > 0)
    % It has both negative and positive
    % Find ratio of negative to positive
    ratio = abs(lims(1)) / (abs(lims(1)) + lims(2));
    neglen = round(m*ratio);
    poslen = m - neglen;
    
    % Just negative
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, neglen);
    newmap1 = zeros(neglen, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap1(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    % Just positive
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, poslen);
    newmap = zeros(poslen, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    % And put 'em together
    newmap = [newmap1; newmap];
    
elseif lims(1) >= 0
    % Just positive
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
else
    % Just negative
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
end