function plotOnSphere(x,y,z,varargin)

    %// Vectors representing each point
    xyz = [x(:), y(:), z(:)].';  %'

    %// One vector of the "first" points and one of the "next" points
    v1 = xyz(:, 1:end-1);
    v2 = xyz(:, 2:end);

    %// Cross product between the vectors of one point and the next
    cv1v2 = cross(v1, v2);

    %// Compute unit vector in the plane defined by v1 and v2
    cv1v2v1 = cross(cv1v2, v1);
    for ii = 1:size(cv1v2v1,2)
        v3(:,ii) = cv1v2v1(:,ii)./norm(cv1v2v1(:,ii));
    end

    %// Figure out the range of the inner angle between v1 and v2
    nc = sqrt(sum(cv1v2.^2, 1));
    t = atan2(nc, dot(v1, v2, 1));

    %// Number of points to sample between any two points on the sphere
    nPoints = 100;

    %// Compute the interpolant
    V = zeros([nPoints, fliplr(size(v1))]);
    for k = 1:numel(t)
        T = linspace(0, t(k), 100);
        V(:,k,:) = (v1(:,k) * cos(T) + v3(:,k) * sin(T)).';
    end

    %// Break the result out into x,y,z parts
    xx = V(:,:,1);
    yy = V(:,:,2);
    zz = V(:,:,3);

    %// Plot the lines
    h = plot3(xx(:), yy(:), zz(:), '-', varargin{:});
    hold on

end


    