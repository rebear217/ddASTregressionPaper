function DZ = D(Z)

    n = length(Z);
    h = 1/(n-1);
    %diffusion in 1-d with Neumann conditions
    %Q = zeros(size(Z));
    %Q(1) = 2*(Z(2)-Z(1));
    %Q(2:end-1) = Z(1:end-2)-2*Z(2:end-1)+Z(3:end);
    %Q(end) = 2*(Z(end-1)-Z(end));
    %Q = Q/h^2;
    %DZ = Q;

    %radial diffusion in 2-d with Neumann conditions
    Qrr = zeros(size(Z));
    Qr = zeros(size(Z));
    r = ((n-1):-1:0)/(n-1);
    r = r';
    Qrr(1) = 2*(Z(2)-Z(1));
    Qrr(2:end-1) = Z(1:end-2)-2*Z(2:end-1)+Z(3:end);
    Qrr(end) = 2*(Z(end-1)-Z(end));
    
    Qr(1) = 0;
    Qr(2:end-1) = (Z(1:end-2)-Z(3:end))/2;
    Qr(end) = 0;

    DZ = Qrr/h^2 + [Qr(1:end-1)./r(1:end-1)/h ; 0];    

end