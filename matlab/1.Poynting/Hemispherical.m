function [Abs_p, Abs_s] = Hemispherical(Abs_TM, Abs_TE, theta_i, d_theta)
% Abs_avg = Abs_TM*0.5 + Abs_TE*0.5;
for index = 1:length(theta_i)
    Abs_temp_p(index) = Abs_TM(index)*cos(theta_i(index)*pi/180)*sin(theta_i(index)*pi/180)*(d_theta*pi/180);
    Abs_temp_s(index) = Abs_TE(index)*cos(theta_i(index)*pi/180)*sin(theta_i(index)*pi/180)*(d_theta*pi/180);
end
Abs_p = 2*sum(Abs_temp_p(1:length(theta_i)));
Abs_s = 2*sum(Abs_temp_s(1:length(theta_i)));
if (Abs_p>(1+1e-6))||(Abs_s>(1+1e-6))
    error('Absorptivity is greater than 1');
end

end
