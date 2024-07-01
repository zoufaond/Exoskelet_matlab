function ref = reference_trajectory(t,t_final,x0,scale)
    zero_to_pi = pi*t/t_final;
    Abd = (sin(zero_to_pi-pi/2)+1)*scale;
    ref = [x0(1:7)+Abd*0;x0(8)+Abd;x0(9:10)+Abd*0];
end