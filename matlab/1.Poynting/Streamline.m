% Richard Z. Zhang 4/30/2014
% Simple vector tracing function from all column vectors

function [x_norm] = Streamline(S_x, S_z, z_norm)

x_norm(1) = 0;
    for ind = 2:length(z_norm)
        x_norm(ind) = x_norm(ind-1) + S_x(:,(ind-1))/S_z(:,(ind-1))*(z_norm(ind)-z_norm(ind-1));
    end
end