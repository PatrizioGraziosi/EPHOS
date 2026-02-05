% this script, as well as the bxsf_to..., must be copied, together the .sh
% file, in the pseudo_inter and pseudo_intra dirs.


modes_info

%
alat = 0.1;
reduce_flag = 'n'; % 'y' if SOC and the system is not magnetic
interpolation_factor = 0;

E_band_window = 0.5; % in eV, range beyong the band edge where the bands are considered, it is related to what you need in transport simulations 

labels_q = {'G', 'Z', 'Y', 'X', 'V', 'U', 'T', 'R'} ;
%


temp = importdata('../freq.dat'); % readmatrix('../freq.dat') *4.136;  % ('../freq.txt') *4.136; % read the frequencies from THz to eV
%temp(isnan(temp)) = 0;
%freq_array = nonzeros(temp);
freq_array = temp;
%disp temp

AA = readmatrix('A_matrix.txt');
a1 = AA(1,:); a2 = AA(2,:); a3 = AA(3,:);


I_c_temp = readmatrix('coord_initial.txt');

for i = size(I_c_temp,1):-1:1
I_c(i,:) = I_c_temp(i,1)*a1 + I_c_temp(i,2) * a2 + I_c_temp(i,3)*a3 ; 
end


shifting_flag = 'no';
[Ek_i_unshifted, ~ ] = bxsf_to_ELECTRA_shifting_option(fileName,material_name,alat,reduce_flag,interpolation_factor,shifting_flag) ;
shifting_flag = 'yes';
[Ek_i_shifted, Fermi]  = bxsf_to_ELECTRA_shifting_option(fileName,material_name,alat,reduce_flag,interpolation_factor,shifting_flag) ;


% definition of the eventual ref. value
Ref_i = max(max(max(Ek_i_unshifted(:,:,:,1)))); % mean(mean(mean(Ek_i_unshifted(:,:,:,1))));
% E_threshold = 0.45; % eV

% definition of VB and CB indexes
[CB_idx_temp, VB_idx_temp, ~ ] = shifting_bands(Ek_i_shifted) ; % the Ek with CB edge set to zero does not interest now

for id_band = (max(size(CB_idx_temp))):-1:1
    E_temp(id_band) = min(min(min(Ek_i_shifted(:,:,:,CB_idx_temp(id_band)))));
end
CB_edge = min(E_temp); 
clear E_temp;
CB_idx = 0*CB_idx_temp;
for id_band = (max(size(CB_idx_temp))):-1:1
    if min(min(min(Ek_i_shifted(:,:,:,CB_idx_temp(id_band))))) - CB_edge < E_band_window
        CB_idx(id_band) = CB_idx_temp(id_band);
    end
end
CB_idx = nonzeros(CB_idx);

for id_band = (max(size(VB_idx_temp))):-1:1
    E_temp(id_band) = max(max(max(Ek_i_shifted(:,:,:,VB_idx_temp(id_band)))));
end
VB_edge = max(E_temp); 
clear E_temp;
VB_idx = 0*VB_idx_temp;
for id_band = (max(size(VB_idx_temp))):-1:1
    if abs(max(max(max(Ek_i_shifted(:,:,:,VB_idx_temp(id_band))))) - VB_edge) < E_band_window
        VB_idx(id_band) = VB_idx_temp(id_band);
    end
end
VB_idx = nonzeros(VB_idx);

if CB_edge > VB_edge
    semiconductor = 'yes';
    Ek_i_shifted_CB = Ek_i_shifted - CB_edge;
    Ek_i_shifted_VB = Ek_i_shifted - VB_edge;
else
    Ek_i_shifted_CB = Ek_i_shifted;
end

% find minima values and positions

% initialization
CB_i = struct();
for id_band = (max(size(CB_idx))):-1:1
    CB_i(id_band).E_ave_unshifted = 0;
    CB_i(id_band).tw = 0;
end
% evaluation
for id_band = max(size(CB_idx)):-1:1

    Ek_temp = Ek_i_shifted_CB(:,:,:,CB_idx(id_band)) ;
    Ek_temp_unshifted = Ek_i_unshifted(:,:,:,CB_idx(id_band)) ;
        
    CB_i(id_band).E_ave_unshifted = mean(mean(mean(Ek_temp_unshifted)));
    CB_i(id_band).tw = 1/4 * ...
        abs(max(max(max(Ek_temp)))-min(min(min(Ek_temp)))); % 2tw*cos...

end

if strcmp(semiconductor,'yes')
    %initialization
    VB_i = struct();
    for id_band = (max(size(VB_idx))):-1:1
        VB_i(id_band).tw = 0;
        VB_i(id_band).E_ave_unshifted_cut = 0;
    end
    % evaluation
    for id_band = max(size(VB_idx)):-1:1
    
        Ek_temp = Ek_i_shifted_VB(:,:,:,VB_idx(id_band)) ;
        Ek_temp_unshifted = Ek_i_unshifted(:,:,:,VB_idx(id_band)) ;
        
        VB_i(id_band).E_ave_unshifted = mean(mean(mean(Ek_temp_unshifted)));
        VB_i(id_band).tw = 1/4 * ...
            abs(max(max(max(Ek_temp)))-min(min(min(Ek_temp)))); % 2tw*cos...
    
    end
end

% initialization of the DP structures
DP_tw = struct();
DP_unshifted = struct();


for im = 1:n_modes
    for iq = 1:n_q
        DP_tw(im,iq).CB = zeros(max(size(CB_idx)), max(size(CB_idx)));
        DP_unshifted(im,iq).CB = zeros(max(size(CB_idx)), max(size(CB_idx)));
        DP_tw(im,iq).freq = 0;
        DP_unshifted(im,iq).freq = 0;
        if strcmp(semiconductor,'yes')
            DP_tw(im,iq).VB = zeros(max(size(VB_idx)), max(size(VB_idx)));
            DP_unshifted(im,iq).VB = zeros(max(size(VB_idx)), max(size(VB_idx)));
        end
        DP_tw(im,iq).d0 = 0;
        DP_unshifted(im,iq).d0 = 0;
    end
end

counter = 0;
for iq = 1:n_q
    for im = 1:n_modes
        counter = counter + 1;
        if counter <= 9
             path_to_dir = ['mod-00',num2str(counter),'/'];
        elseif counter <= 99
            path_to_dir = ['mod-0',num2str(counter),'/'];
        else
            path_to_dir = ['mod-',num2str(counter),'/'];
        end
             
        cd(path_to_dir) % the bxsf_to... files must be copied in every dir.
        shifting_flag = 'no';
        [Ek_disp_unshifted, ~ ] = bxsf_to_ELECTRA_shifting_option(fileName,material_name,alat,reduce_flag,interpolation_factor,shifting_flag) ;
        shifting_flag = 'yes';
        [Ek_disp_shifted, Fermi]  = bxsf_to_ELECTRA_shifting_option(fileName,material_name,alat,reduce_flag,interpolation_factor,shifting_flag) ; %#ok<*ASGLU>
        
        Ref_disp = max(max(max(Ek_disp_unshifted(:,:,:,1)))); % mean(mean(mean(Ek_disp_unshifted(:,:,:,1))));


        M_c_temp = readmatrix('coord.txt'); % read the modulated coord
        for i = size(M_c_temp,1):-1:1
            M_c(i,:) = M_c_temp(i,1)*a1 + M_c_temp(i,2) * a2 + M_c_temp(i,3)*a3 ; 
        end
        
%        D_c = M_c - I_c;
%        displ = sqrt( D_c(:,1).^2 + D_c(:,2).^2 + D_c(:,3).^2 );
	displ =  sqrt( sum( (M_c-I_c).^2,2) );
        N = displ > 10*mean(displ);
disp(displ(N)')
	displ(N) = max(max(AA)) - displ(N) ; % 1-displ(N);
	d0 = mean(displ)  % average displacement in Angstrom


        % initialize the structure(s)
        CB_disp = struct();
        for id_band = (max(size(CB_idx))):-1:1
            CB_disp(id_band).E_ave_unshifted = 0;
            CB_disp(id_band).tw = 0;
        end

        for id_band = (max(size(CB_idx_temp))):-1:1
            E_temp(id_band) = min(min(min(Ek_disp_shifted(:,:,:,CB_idx_temp(id_band)))));
        end
        CB_edge_temp = min(E_temp);
        clear E_temp;
        for id_band = (max(size(VB_idx_temp))):-1:1
            E_temp(id_band) = max(max(max(Ek_disp_shifted(:,:,:,VB_idx_temp(id_band)))));
        end
        VB_edge_temp = max(E_temp);
        if strcmp(semiconductor,'yes')
            Ek_disp_shifted_CB = Ek_disp_shifted - CB_edge_temp;
            Ek_disp_shifted_VB = Ek_disp_shifted - VB_edge_temp;

            VB_disp = struct();
            for id_band = (max(size(VB_idx))):-1:1
                VB_disp(id_band).E_ave_unshifted = 0;
                VB_disp(id_band).tw = 0;
            end

        else
            Ek_disp_shifted_CB = Ek_disp_shifted;
        end


        % the part below, until line 214, must be repeated for VB
        for id_band = max(size(CB_idx)):-1:1
            Ek_temp = Ek_disp_shifted_CB(:,:,:,CB_idx(id_band)) ;
            Ek_temp_unshifted = Ek_disp_unshifted(:,:,:,CB_idx(id_band)) ;
            
            
            CB_disp(id_band).E_ave_unshifted =  mean(mean(mean(Ek_temp_unshifted)));
            CB_disp(id_band).tw = 1/4 * ...
                abs(max(max(max(Ek_temp)))-min(min(min(Ek_temp)))); 

        end


        % compute the deformation potentials DP
        for id_band = max(size(CB_idx)):-1:1

            DP_tw(im,iq).CB(id_band,id_band) = ...
                abs( CB_disp(id_band).tw - CB_i(id_band).tw) / d0 ; % matrixes nb x nb

            DP_unshifted(im,iq).CB(id_band,id_band) = ...
                abs( CB_disp(id_band).E_ave_unshifted - Ref_disp - ...
                ( CB_i(id_band).E_ave_unshifted - Ref_i ) ) / d0 ; 


            for id_band_2 = max(size(CB_idx)):-1:1 % inter-band evaluated from Davydov split

                if ne(id_band,id_band_2)

                    DP_tw(im,iq).CB(id_band,id_band_2) = ...
                    abs( abs( CB_disp(id_band).E_ave_unshifted - CB_disp(id_band_2).E_ave_unshifted ) ...
                    - abs( CB_i(id_band).E_ave_unshifted - CB_i(id_band_2).E_ave_unshifted ) )/ d0 ;
                    DP_tw(im,iq).CB(id_band_2,id_band) = DP_tw(im,iq).CB(id_band,id_band_2);
    
                    DP_unshifted(im,iq).CB(id_band,id_band_2) = DP_tw(im,iq).CB(id_band,id_band_2);
                    DP_unshifted(im,iq).CB(id_band_2,id_band) = DP_unshifted(im,iq).CB(id_band,id_band_2);

                end

            end

        end


        % as above for VB, from line 178
        if strcmp(semiconductor,'yes')

            for id_band = max(size(VB_idx)):-1:1

                Ek_temp = Ek_disp_shifted_VB(:,:,:,VB_idx(id_band)) ;
                Ek_temp_unshifted = Ek_disp_unshifted(:,:,:,VB_idx(id_band)) ;

                VB_disp(id_band).E_ave_unshifted =  mean(mean(mean(Ek_temp_unshifted)));
                VB_disp(id_band).tw = 1/4 * ...
                abs(max(max(max(Ek_temp)))-min(min(min(Ek_temp)))); 

             end


            % compute the deformation potentials DP
            for id_band = max(size(VB_idx)):-1:1
    
                DP_tw(im,iq).VB(id_band,id_band) = ...
                    abs( VB_disp(id_band).tw - VB_i(id_band).tw) / d0 ; % matrixes nb x nb
    
                DP_unshifted(im,iq).CB(id_band,id_band) = ...
                    abs( VB_disp(id_band).E_ave_unshifted - Ref_disp - ...
                    ( VB_i(id_band).E_ave_unshifted - Ref_i ) ) / d0 ; 
    
    
                for id_band_2 = max(size(VB_idx)):-1:1 % inter-band evaluated from Davydov split
                    if ne(id_band,id_band_2)

                        DP_tw(im,iq).VB(id_band,id_band_2) = ...
                        abs( abs( VB_disp(id_band).E_ave_unshifted - VB_disp(id_band_2).E_ave_unshifted ) ...
                        - abs( VB_i(id_band).E_ave_unshifted - VB_i(id_band_2).E_ave_unshifted ) )/ d0 ;
                        DP_tw(im,iq).VB(id_band_2,id_band) = DP_tw(im,iq).VB(id_band,id_band_2);
        
                        DP_unshifted(im,iq).VB(id_band,id_band_2) = DP_tw(im,iq).VB(id_band,id_band_2);
                        DP_unshifted(im,iq).VB(id_band_2,id_band) = DP_unshifted(im,iq).VB(id_band,id_band_2);

                    end
    
                end

            end        
        end


        save mode_point_specific_data
        cd ../
        
%        whos freq_array
        f_counter = (iq-1)*total_modes + im + modes_to_skip ;
        DP_tw(im,iq).freq = freq_array(f_counter);
        DP_unshifted(im,iq).freq = freq_array(f_counter);
        DP_tw(im,iq).d0 = d0;
        DP_unshifted(im,iq).d0 = d0;

    end
end



% for the mod-all case
cd mod-all/ % the bxsf_to... files must be copied in every dir.
shifting_flag = 'no';
[Ek_disp_unshifted, ~ ] = bxsf_to_ELECTRA_shifting_option(fileName,material_name,alat,reduce_flag,interpolation_factor,shifting_flag) ;
shifting_flag = 'yes';
[Ek_disp_shifted, Fermi]  = bxsf_to_ELECTRA_shifting_option(fileName,material_name,alat,reduce_flag,interpolation_factor,shifting_flag) ;

Ref_disp = max(max(max(Ek_disp_unshifted(:,:,:,1)))); % mean(mean(mean(Ek_disp_unshifted(:,:,:,1))));

M_c = readmatrix('coord.txt'); % read the modulated coord
%D_c = M_c - I_c;
%displ = sqrt( D_c(:,1).^2 + D_c(:,2).^2 + D_c(:,3).^2 );
displ =  sqrt( sum( (M_c-I_c).^2,2) );
N = displ > 10*mean(displ);
disp(displ(N)')
displ(N) = max(max(AA)) - displ(N) ; % 1-displ(N);
d0 = mean(displ) ; % average displacement in Angstrom


% initialize the structure(s)
CB_disp = struct();
for id_band = (max(size(CB_idx))):-1:1
    CB_disp(id_band).E_ave_unshifted = 0;
    CB_disp(id_band).tw = 0;
end

for id_band = (max(size(CB_idx_temp))):-1:1
    E_temp(id_band) = min(min(min(Ek_disp_shifted(:,:,:,CB_idx_temp(id_band)))));
end
CB_edge_temp = min(E_temp);
clear E_temp;
for id_band = (max(size(VB_idx_temp))):-1:1
    E_temp(id_band) = max(max(max(Ek_disp_shifted(:,:,:,VB_idx_temp(id_band)))));
end
VB_edge_temp = max(E_temp);
if strcmp(semiconductor,'yes')
    Ek_disp_shifted_CB = Ek_disp_shifted - CB_edge_temp;
    Ek_disp_shifted_VB = Ek_disp_shifted - VB_edge_temp;

    VB_disp = struct();
    for id_band = (max(size(VB_idx))):-1:1
        VB_disp(id_band).E_ave_unshifted = 0;
        VB_disp(id_band).tw = 0;
    end

else
    Ek_disp_shifted_CB = Ek_disp_shifted;
end


% the part below, until line 214, must be repeated for VB
for id_band = max(size(CB_idx)):-1:1
    Ek_temp = Ek_disp_shifted_CB(:,:,:,CB_idx(id_band)) ;
    Ek_temp_unshifted = Ek_disp_unshifted(:,:,:,CB_idx(id_band)) ;
    
    
    CB_disp(id_band).E_ave_unshifted =  mean(mean(mean(Ek_temp_unshifted)));
    CB_disp(id_band).tw = 1/4 * ...
        abs(max(max(max(Ek_temp)))-min(min(min(Ek_temp)))); 

end


% compute the deformation potentials DP
for id_band = max(size(CB_idx)):-1:1

    all_DP_tw.CB(id_band,id_band) = ...
        abs( CB_disp(id_band).tw - CB_i(id_band).tw) / d0 ; % matrixes nb x nb

    all_DP_unshifted.CB(id_band,id_band) = ...
        abs( CB_disp(id_band).E_ave_unshifted - Ref_disp - ...
        ( CB_i(id_band).E_ave_unshifted - Ref_i ) ) / d0 ; 


    for id_band_2 = max(size(CB_idx)):-1:1 % inter-band evaluated from Davydov split

        if ne(id_band,id_band_2)

            all_DP_tw.CB(id_band,id_band_2) = ...
            abs( abs( CB_disp(id_band).E_ave_unshifted - CB_disp(id_band_2).E_ave_unshifted ) ...
            - abs( CB_i(id_band).E_ave_unshifted - CB_i(id_band_2).E_ave_unshifted ) )/ d0 ;
             all_DP_tw.CB(id_band_2,id_band) = all_DP_tw.CB(id_band,id_band_2);

            all_DP_unshifted.CB(id_band,id_band_2) = all_DP_tw.CB(id_band,id_band_2);
            all_DP_unshifted.CB(id_band_2,id_band) = all_DP_tw.CB(id_band,id_band_2);

        end

    end

end


% as above for VB, from line 178
if strcmp(semiconductor,'yes')

    for id_band = max(size(VB_idx)):-1:1

        Ek_temp = Ek_disp_shifted_VB(:,:,:,VB_idx(id_band)) ;
        Ek_temp_unshifted = Ek_disp_unshifted(:,:,:,VB_idx(id_band)) ;

        VB_disp(id_band).E_ave_unshifted =  mean(mean(mean(Ek_temp_unshifted)));
        VB_disp(id_band).tw = 1/4 * ...
        abs(max(max(max(Ek_temp)))-min(min(min(Ek_temp)))); 

     end


    % compute the deformation potentials DP
    for id_band = max(size(VB_idx)):-1:1

         all_DP_tw.VB(id_band,id_band) = ...
            abs( VB_disp(id_band).tw - VB_i(id_band).tw) / d0 ; % matrixes nb x nb

         all_DP_unshifted.VB(id_band,id_band) = ...
            abs( VB_disp(id_band).E_ave_unshifted - Ref_disp - ...
            ( VB_i(id_band).E_ave_unshifted - Ref_i ) ) / d0 ; 


        for id_band_2 = max(size(VB_idx)):-1:1 % inter-band evaluated from Davydov split
            if ne(id_band,id_band_2)

                all_DP_tw.VB(id_band,id_band_2) = ...
                abs( abs( VB_disp(id_band).E_ave_unshifted - VB_disp(id_band_2).E_ave_unshifted ) ...
                - abs( VB_i(id_band).E_ave_unshifted - VB_i(id_band_2).E_ave_unshifted ) )/ d0 ;
                all_DP_tw.VB(id_band_2,id_band) = all_DP_tw.VB(id_band,id_band_2);

                all_DP_unshifted.VB(id_band,id_band_2) = all_DP_tw.VB(id_band,id_band_2);
                all_DP_unshifted.VB(id_band_2,id_band) = all_DP_unshifted.VB(id_band,id_band_2);

            end

        end

    end
    all_DP_tw.d0 = d0;
    all_DP_unshifted.d0 = d0;

end
cd ..

save('DeformationPotentials_data', "DP_tw", "DP_unshifted", "all_DP_tw", "all_DP_unshifted")


n_CB = max(size(CB_idx));
n_CB_inter = (n_CB^2-n_CB)/2;
n_VB = max(size(VB_idx));
n_VB_inter = (n_VB^2-n_VB)/2;
for i_b = 1:n_CB % intra-band CB
    for im = n_modes:-1:1
        for iq = n_q:-1:1         
            matrix_f(im,iq) = DP_tw(im,iq).freq ;
            matrix(im,iq) = DP_tw(im,iq).CB(i_b,i_b) ;
        end
    end
    table_csv = array2table( matrix_f, 'VariableNames',  labels_q ) ;
    csv_name = 'freqs.csv';
    writetable( table_csv, csv_name ) ;
    table_csv = array2table( matrix, 'VariableNames',  labels_q ) ;
    csv_name = ['tw_mod_CB',i_b,'.csv'] ;
    writetable( table_csv, csv_name ) ;
end
for i_b = 1:n_VB % intra-band VB
    for im = n_modes:-1:1
        for iq = n_q:-1:1         
            matrix(im,iq) = DP_tw(im,iq).VB(i_b,i_b) ;
        end
    end
    table_csv = array2table( matrix, 'VariableNames',  labels_q ) ;
    csv_name = ['tw_mod_VB',i_b,'.csv'] ;
    writetable( table_csv, csv_name ) ;
end
for i_b = 1:n_CB % inter-band CB
    for i_b2 = 1:n_CB
        if ne (i_b,i_b)
            for im = n_modes:-1:1
                for iq = n_q:-1:1         
                    matrix(im,iq) = DP_tw(im,iq).CB(i_b,i_b2) ;
                end
            end
            table_csv = array2table( matrix, 'VariableNames',  labels_q ) ;
            csv_name = ['DavSplit_mod_CB',i_b,'_',i_b2,'.csv'] ;
            writetable( table_csv, csv_name ) ;
        end
    end
end
for i_b = 1:n_VB % inter-band VB
    for i_b2 = 1:n_VB
        if ne (i_b,i_b)
            for im = n_modes:-1:1
                for iq = n_q:-1:1         
                    matrix(im,iq) = DP_tw(im,iq).VB(i_b,i_b2) ;
                end
            end
            table_csv = array2table( matrix, 'VariableNames',  labels_q ) ;
            csv_name = ['DavSplit_mod_VB',i_b,'_',i_b2,'.csv'] ;
            writetable( table_csv, csv_name ) ;
        end
    end
end




% % subfunction
% function [pseudoconduction_bands, pseudovalence_bands, Ek] = shifting_bands(Ek) % %#codegen
% % shifting of the bands to zero and identify the CB and VB indexes
% 
%     Emaximalextreme = zeros(1,size(Ek,4)) ;
%     Eminimalextreme = zeros(1,size(Ek,4)) ;
%     for id_n = 1:size(Ek,4)
%         Emaximalextreme(id_n) = max(max(max(Ek(:,:,:,id_n))));
%         Eminimalextreme(id_n) = min(min(min(Ek(:,:,:,id_n))));
%     end
% 
%     [~,pseudoconduction_bands] = find( abs(Emaximalextreme) > abs(Eminimalextreme) );
%     [~,pseudovalence_bands] = find( abs(Emaximalextreme) < abs(Eminimalextreme) );
%     
%     if numel(pseudoconduction_bands)>0
%         new_Fermi = min(min(min(min(Ek(:,:,:,pseudoconduction_bands(1):pseudoconduction_bands(size(pseudoconduction_bands,2)))))));
%         Ek = Ek - new_Fermi; % this the conduction band, or the flipped valence band, start from zero, positive values of the EF_array are into the band and negative EF values are into the gap
%     end
% end
