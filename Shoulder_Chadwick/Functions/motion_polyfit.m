function out = motion_polyfit(motion_full,plot_polyfits)
time = motion_full(:,1)-motion_full(1,1);
mot_eul = motion_full(:,8:end-1)*pi/180;
out.tout = time;
for jnt = 1:length(mot_eul(1,:))
    p = polyfit(time,mot_eul(:,jnt),6);
    polyvals(:,jnt) = p;
    out.data(:,jnt) = polyval(p,time);
    if strcmp(plot_polyfits,'plot_polyfits')
        figure
        plot(time,mot_eul(:,jnt),'o')
        hold on
        plot(time,out.data(:,jnt))
    end
end

end