function [RingingPeriod] = getRingingPeriod(SamplingRatio)
%GETRINGINGPERIOD Summary of this function goes here
%   Detailed explanation goes here
    RingingPeriod = 2/sqrt(SamplingRatio); % Unit: pixel
end

