function [percent_of_flies_to_respond] = calculate_percent(path,th_angle)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% Read the CSV file into a table
data = readtable(path); 
data = table2array(data(:,2:end));
if th_angle ~= false
    numerator = (isnan(data) == false) & (data ~= 999) & (data ~= 1000) & (abs(data) < th_angle);
else
    numerator = (isnan(data) == false) & (data ~= 999) & (data ~= 1000);
end
denominator = (isnan(data) == false) & (data ~= 999);
percent_of_flies_to_respond = sum(numerator, 1) * 100 ./ sum(denominator,1);


end