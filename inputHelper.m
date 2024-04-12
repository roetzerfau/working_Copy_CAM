% g = @(x) 2000/1061.8 * (x - 365*cos(2*pi*x/365)/(2*pi) + 365/(2*pi));

% g = @(x) 2000/765.46 * (3/4*x - 365/(8*pi)*cos(2*pi*x/365) + 365/(8*pi));
g = @(x) 500/1030.9 * (x - 365/(4*pi)*cos(2*pi*x/365) + 365/(4*pi));

inputVectorHelper = zeros(1000,1);

for i = 1:1000
    inputVectorHelper(i) = floor(g(i));
end

inputVector = zeros(1000,1);

inputVector(1) = inputVectorHelper(1);
inputVector(2:end) = inputVectorHelper(2:end) - inputVectorHelper(1:end-1);
