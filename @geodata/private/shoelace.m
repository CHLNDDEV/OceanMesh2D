function area = shoelace(x,y)
n = length(x);
xp = [x; x(1)];
yp = [y; y(1)];
area = 0;
for i = 1:n
    area = area + det([xp(i), xp(i+1); yp(i), yp(i+1)]);
end
area = 1/2*abs(area);