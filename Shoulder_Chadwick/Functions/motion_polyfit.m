function [time_new, data] = motion_polyfit(motion_full,plot_polyfits,numdata)
time = motion_full(:,1)-motion_full(1,1);
time_new = linspace(0,time(end,1),numdata);
mot_eul = motion_full(:,8:end)*pi/180;
for jnt = 1:length(mot_eul(1,:))
    p = polyfit(time,mot_eul(:,jnt),4);
    polyvals(:,jnt) = p;
    data(:,jnt) = polyval(p,time_new);
    if strcmp(plot_polyfits,'plot_polyfits')
        figure
        plot(time,mot_eul(:,jnt),'o')
        hold on
        plot(time_new,data(:,jnt))
    end
end

end