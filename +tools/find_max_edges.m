function [edges,theta_max] = find_max_edges(array,options)
    arguments
        array (:,:) double {mustBeNumeric}
        options.res (1,1) int16 {mustBePositive} = 100
    end
    
    Gx = [[-1,0,1];[-2,0,2];[-1,0,1]];
    Gy = [[1,2,1];[0,0,0];[-1,-2,-1]];
    
    conv_x = abs(conv2(Gx,array));
    conv_y = abs(conv2(Gy,array));
    conv_x = conv_x(3:end-2,3:end-2);
    conv_y = conv_y(3:end-2,3:end-2);
    
    dtheta = 2*pi/double(options.res);
    thetas = -pi:dtheta:pi;
    sums = zeros(size(thetas));
    for itheta = 1:length(thetas)
        theta = thetas(itheta);
        sums(itheta) = sum(cos(theta)*conv_x + sin(theta)*conv_y,'all');
    end
    [~,imax] = max(sums);
    theta_max = thetas(imax);
    edges = cos(theta_max)*conv_x + sin(theta_max)*conv_y;
end