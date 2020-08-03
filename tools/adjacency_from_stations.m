function [A,D] = adjacency_from_stations(stations_data, max_km)
% Constants
COOR_TO_MILES = 60;
MILES_TO_KM = 1.852;

lat_long = [stations_data{3} stations_data{4}];
A = zeros(length(stations_data{1}));
D = zeros(length(stations_data{1}));
for i=1:length(A)
    pos1 = [lat_long(i,1) lat_long(i,2)];
    for j=i:length(A)
        pos2 = [lat_long(j,1) lat_long(j,2)];
        D(j,i) = distance(pos1,pos2)*COOR_TO_MILES*MILES_TO_KM;
        if D(j,i) <= max_km && i ~= j
            A(j,i) = 1;
        else
            A(j,i) = 0;
        end
    end
end
A = A+A';