motion_coords = readmatrix('motion.txt');
motion_coords = motion_coords(:,1:10);
motion_coords(:,3) = motion_coords(:,3)*(-1); % kvuli zmene poradi rotace
motion_coords(:,7) = motion_coords(:,7)*(-1); % kvuli zmene poradi rotace
motion_coords(:,8) = motion_coords(:,9)*(-1); % v opensimu je zde rotace zyz, pouze jsem elevaci prepsal na rotaci kolem X v simscapu
motion_coords(:,9) = motion_coords(:,9)*(0); %zbaveni druhych dvou rotaci u humeru
motion_coords(:,10) = motion_coords(:,10)*(0); %zbaveni druhych dvou rotaci u humeru

time = linspace(0,3.88,195);
p_coords = zeros(9,5);
for i=2:10
    p_coords(i-1,:) = polyfit(time,motion_coords(:,i),4);
end
