clear; 
close all; 
clc;
% Define variables
number_cracks = 100;
divisions = 35;
length_whole_segment = 1;
gauss_per_division = 2;

%% This is a new addition
x = x + 1


% Space is a vector of the endpoints of all the segments
space = 0:length_whole_segment/divisions:length_whole_segment;
numberGauss = gauss_per_division*divisions;
numbered_divisions = 1:divisions;

% Create connectivity matrix
connect = [1:divisions-1;2:divisions];

% Caclulate length of each segment, using connectivity matrix
length_each_segment = space(connect(2,:)) - space(connect(1,:));

random_cracks = rand(number_cracks,1);

cracks_in_division = zeros(divisions,number_cracks);
%this comment is by Farah
%looks like one for loop in the code. put this loop inside a function. so
%the input of that function is
%divisions,number_cracks,random_cracks,space and output is
%cracks_in_divisions
for i=1:divisions
    for j=1:number_cracks
        if random_cracks(j)>=space(i) && random_cracks(j)<space(i+1)
            cracks_in_division(i,j) = random_cracks(j);
            j = j+1;
        else
            j = j+1;
        end
    end
    i = i+1;
end

% Show the number of cracks per division
%if sum(cracks_each_division) == number_cracks
%   disp('All cracks counted')
%else
%     disp('Cracks omitted/double counted')
% end


for i=1:divisions
    for j=1:number_cracks
        [r,c,v] = find(cracks_in_division(i,:));
        if isempty(v)
            
        else
            for k=1:length(v)
                s1(i,k) = (v(k)-space(i+1))./(space(i)-space(i+1));
                s2(i,k) = (v(k)-space(i))./(space(i+1)-space(i));
            end
        end
        
    end
    i = i+1;
end

% Check to make sure s1+s2 is 1 (or zero, for empty values)
check_matrix = s1+s2;


s1_mod = zeros(divisions+1,1);
s2_mod = zeros(divisions+1,1);

s1_mod(1:divisions,1) = sum(s1,2);
s2_mod(2:divisions+1,1) = sum(s2,2);

% Calculate and output the nodal weight
nodal_weight = (s1_mod+s2_mod)/number_cracks

% Check to make sure all cracks are accounted for
% Note that there may be a slight deviation from 1, just simply because of
% calculation errors in MATLAB, stemming from the random values.
% If the value is within .00000001, I feel confident that the calculation
% is correct...
if sum(nodal_weight)-1 < 1e-8
    disp('All points accounted for :)')
end

%%% This is totally wrong!!!!!!!!!!!!!!!!
%%%======================================
stdgw(1)=1/sqrt(3);
stdgw(2)=-1/sqrt(3);
S1_g=rg-1;
S2_g=1-S1_g;









