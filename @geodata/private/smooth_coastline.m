function [segment,closed] = smooth_coastline(segment,window,plot_on) 
% smooth polygons and coastline by applying window pt moving average 
% kjr apr. 2017
% modified by WJP June 23 2017
Isnan = find(isnan(segment(:,1)));
Isnan = vertcat(0,Isnan); 
closed = 1;
if length(Isnan) - 1 == 0
    iseg = segment;
    if iseg(1,1) == iseg(end,1)
        % Is polygon so need pad to make cyclic
        iseg = [iseg(end-floor(window/2):end,:); iseg; ...
                iseg(1:floor(window/2),:)];
        if length(iseg) > window
          %iseg(:,1) = smooth(iseg(:,1),window);
          iseg(:,1)= fastsmooth(iseg(:,1),window,1,1); 
          %iseg(:,2) = smooth(iseg(:,2),window);
          iseg(:,2)= fastsmooth(iseg(:,2),window,1,1); 
        end
        segment = iseg(floor(window/2)+1:end-floor(window/2)-1,:);
        segment(end,:) = segment(1,:);
    else
        closed = 0;
        if length(iseg) > window 
          %segment(:,1) = smooth(segment(:,1),window);
           segment(:,1)= fastsmooth(segment(:,1),window,1,1); 
          %segment(:,2) = smooth(segment(:,2),window);
           segment(:,2)= fastsmooth(segment(:,1),window,1,1); 
        end
    end
else
    for i = 1 : length(Isnan) - 1
        iseg = segment(Isnan(i)+1:Isnan(i+1)-1,:);
        if ~isempty(iseg)
            if iseg(1,1) == iseg(end,1) && iseg(1,2) == iseg(end,2)
                % Is polygon so need pad to make cyclic
                iseg = [iseg(end-floor(window/2):end,:); iseg; ...
                        iseg(1:floor(window/2),:)];
                if length(iseg) > window
                  %iseg(:,1) = smooth(iseg(:,1),window);
                  iseg(:,1)= fastsmooth(iseg(:,1),window,1,1); 
                  %iseg(:,2) = smooth(iseg(:,2),window);
                  iseg(:,2)= fastsmooth(iseg(:,2),window,1,1); 
                end
                segment(Isnan(i)+1:Isnan(i+1)-1,:) = ...
                       iseg(floor(window/2)+1:end-floor(window/2)-1,:);
                segment(Isnan(i+1)-1,:) = segment(Isnan(i)+1,:);
            else
                closed = 0;
                if length(iseg) > window
                  %segment(Isnan(i)+1:Isnan(i+1)-1,1) = smooth(iseg(:,1),window);
                  segment(Isnan(i)+1:Isnan(i+1)-1,1) = fastsmooth(iseg(:,1),window,1,1);
                  %segment(Isnan(i)+1:Isnan(i+1)-1,2) = smooth(iseg(:,2),window); 
                  segment(Isnan(i)+1:Isnan(i+1)-1,2) = fastsmooth(iseg(:,2),window,1,1); 
                end

            end
        end
    end
end
% Plot the smoothed coastline
if plot_on >= 1
    hold on;
    plot(segment(:,1),segment(:,2),'.-');
end

end
