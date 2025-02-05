        % qpoints = {'0 0 0', '-1/2 1/2 0', '0 1/2 0' , '-1/2 0 0',  '0 0 1/2'}; %
% G, A, Z, X, Y  P2_1/c
% qpoints = {'0 0 0', '0 0 1/2', '0 1/2 0' , '1/2 0 0',  '1/2 1/2 0', ...
%     '1/2 0 1/2', '0 1/2 1/2', '1/2 1/2 1/2'}; % G, Z, Y, X, V, U, T, R - bilbao P-1

qpoints = { '0 0 0',  '0 0.5 0',    '0 0 0.5',   '0.5 0 0',  ...
    '0.5 0.5 0',  '0.5 0 0.5',  '0 0.5 0.5',    '0.5 0 0.5',  ...
    '0.5 0.5 0.5'} ; % P2_1, s.g. 14

% qpoints = { '0 0 0', '0 0.5 0', '0 0 0.5',  '0 0.5 0.5', ...
%     '0.5 0 0', '0.5 0.5 0', '-0.5 0 0.5',  '-0.5 0.5 0.5'}; % P2_1, s.g. 4, 
% % G        Z          B           D             Y        C             A             E  
% qpoints = { '0 0 0',  '0 0 0.5'};

% 
% qpoints = { '0 0 0', '0.5 0.5 0', '0.5 0.5 0.5',  '0 0 0.5', ...
%     '0 0.5 0', '0 0.5 0.5'}; % Cmce, s.g. 64, 
% %                G         Y           T              Z             
% %        S        R             

RANDOM = 10; % 0 to take all

max_mode = 40 ; % 23;%  32;    % 66; ROQSAT phBTBT10
initial_mode = 1; %  2;       "



if floor(RANDOM) > 1
    n_r = floor(RANDOM);
    id_r = floor( rand(1,n_r)*size(qpoints,2)*(max_mode-initial_mode+1) );
else
    id_r = 0;
end

final_string = [ ];
random_string = [ ];
c_r = 0;
for i_p = 1:size(qpoints,2)
    for i_m = initial_mode:max_mode        
        final_string = [final_string, ', ',char(qpoints(i_p)), ' ', num2str(i_m), ' ', num2str(0.1)];
        
        c_r = c_r+1;
        if ismember(c_r,id_r)
            random_string = [random_string, ', ',char(qpoints(i_p)), ' ', num2str(i_m), ' ', num2str(0.1)];
        end
    end
end
disp(' ')
disp(random_string)
disp(' ')
disp(final_string)
disp(' ')