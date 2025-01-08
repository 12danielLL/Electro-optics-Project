

function circ_output = circ(L, N, R)
    delta_x=L/N;
    x = linspace(-L/2, L/2 - delta_x, N);
    [X, Y] = meshgrid(x, x); 
    r = sqrt(X.^2 + Y.^2);
    circ_output = r < R;
end
