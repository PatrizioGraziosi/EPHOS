% get a comprehensive Def- Pot. and write scattering matrixes in .mat file
% load('Tetracene_P1_DeformationPotentials_elaborated_01_all_bxsf_same_size.mat')
% save_filename = 'scattering_parameters_Tetra_P1_intra1.mat';

% load('DeformationPotentials_data_P1_inter_002_homog_bxsf.mat')
% save_filename = 'temp.mat';
% load('Tetracene_film_DeformationPotentials_elaborates_Tetracene_film_1bis.mat')
% save_filename = 'scattering_parameters_Tetra_film_1.mat';
clearvars


save_filename = 'scattering_parameters_TIPS_II.mat' ; % 'scattering_parameters_C8BTBTC8_sg2.mat' ;% 'scattering_parameters_phBTBT10_ROQSAT_G100_trials_1.mat'; % 'scattering_parameters_phBTBT10_herringbone_transp_17277_840irr.mat'; % 'scattering_parameters_phBTBT10_ROQSAT_VdW.mat';

EPC_filename = 'TIPS_II_DeformationPotentials_data' ; % 'C8BTBTC8_triclinic_DeformationPotentials_data.mat';% 'ROQSAT_G100_DeformationPotentials_data_modified.mat' ; % 'ROQSAT_G100_DeformationPotentials_data.mat';%'rubrene_ortho_02tag_443mesg_DeformationPotentials_data.mat' ; % 'phBTBT10_ROQSAT_DeformationPotentials_data_1mode_removed.mat' ; % phBTBT10_herringbone_DeformationPotentials_4modes_removed.mat' ; %  'phBTBT10_ROQSAT_DeformationPotentials_data.mat'; %'Anthra_1_DeformationPotentials_data.mat';

phononDOS_filename = 'TIPS_II_total_dos.dat' ; %  'C8BTBTC8_P-1_total_dos.dat'; % 'phBTBT10_herringbone_total_dos.dat' ; % 'phBTBT10_ROQSAT_VdW_total_dos.dat';

% n_m = 42 ; % 65 ;% 28;
% n_q = 6 ; % 9 ;% 8;

f_cut = 2.16; %  5.5; % 2.8 ; % 3.2 ; % 2.7 ; % phBTBT10 herr. 2.7, ROQSAT 3.2 ; % % Anthra4.2 ; % Naptha 4.7; %TIPS %  4 fot HT, 4.5 for LT, 5 for TF

rho_mass_density = 1.06e3; % 1.2255e3 ; % 1.293e3 ; % 1.277e3 ; % 1.66*0.78e3  ; % 0.73e3 ; %  0.62e3 ; % 0.73e3 1.14e3 ; % 1.3e3 ;


also_intra = 'no';
n_m_intra_array = 4; % array of number of grouped intra- phonons, e.g. [4, 2]
f_cut_intra = 5.5;
intra_filename = 'P1_intra_DeformationPotentials_data_elab.mat';

remove_Gamma_ac = 'no';


% hbar_w_eff_inter = 2.25 *4.136e-3 ; % Tetra P1 2.27, Naphtha 2.14, 2.25 film
% 
% hbar_w_eff_intra = 5.05 *4.136e-3 ;


load(EPC_filename)

bands_type_CB =  1:size(DP_tw(1,1).CB,1);
bands_type_VB =  1:size(DP_tw(1,1).VB,1);

n_m = size(DP_tw,1);
n_q = size(DP_tw,2);

n_process_odp_CB = n_m; % these should be n_m + length of n_m_intra_array
n_process_ivs_CB = n_m;

n_process_odp_VB = n_m;
n_process_ivs_VB = n_m;





imported = importdata(phononDOS_filename);
DOS_ph = imported.data;



fr = DOS_ph(DOS_ph(:,1)<f_cut,1);
dosr = DOS_ph(DOS_ph(:,1)<f_cut,2);

n_CB = max(size(DP_tw(1,1).CB)) ;
n_VB = max(size(DP_tw(1,1).VB)) ;

D_odp_e_M = zeros(n_process_odp_CB,n_CB);
hbar_w_odp_e_M = zeros(n_process_odp_CB,n_CB);

D_ivs_e_M = zeros(n_CB,n_CB,n_process_ivs_CB);
hbar_w_ivs_e_M = zeros(n_CB,n_CB,n_process_ivs_CB);

D_odp_h_M = zeros(n_process_odp_VB,n_VB);
hbar_w_odp_h_M = zeros(n_process_odp_VB,n_VB);

D_ivs_h_M = zeros(n_VB,n_VB,n_process_ivs_VB);
hbar_w_ivs_h_M = zeros(n_VB,n_VB,n_process_ivs_VB);



for i_b = 1:n_CB
    
    for im = n_m :-1:1    
        for iq = n_q :-1:1

        f_t = DP_tw(im,iq).freq;

        DOS_t= interp1(fr,dosr,f_t, 'makima');
        if DOS_t < 0 || f_t < 0
            DOS_t = 0;
        end
        DOS_M(im,iq) = DOS_t;

        if strcmp(remove_Gamma_ac,'yes') && im < 4 && iq == 1
            DOS_M(iq,im) = 0;
        end

        DP_temp(im,iq) = DP_tw(im,iq).CB(i_b,i_b) ;
        
        end
    end


    for im = n_m:-1:1
        for iq = n_q:-1:1
            freq_temp(im,iq) = DP_tw(im,iq).freq ;
        end
    end
    w_v2_array = ( (sum(freq_temp.*DOS_M,2))./sum(DOS_M,2) )*4.136e-3 ; % in eV



     % D_odp_e_M(1,i_b) = ...
     %     sqrt(sum(sum(DP_temp.^2.*DOS_M))/(sum(sum(DOS_M)))) *1e10;

     D_odp_e_M(:,i_b) = ...
         sqrt(...
         sum(DP_temp.^2.*DOS_M,2)./(sum(DOS_M,2))...
         ) *1e10;

     hbar_w_odp_e_M(:,i_b) = w_v2_array;
            
     for i_b2 = 1:n_CB
        if ne(i_b,i_b2)
            for im = n_m:-1:1    
                for iq = n_q:-1:1
                    DP_temp(im,iq) =  DP_tw(im,iq).CB(i_b,i_b2) ;
                end
            end

            D_ivs_e_M(i_b,i_b2,:) = ...
                sqrt(...
         sum(DP_temp.^2.*DOS_M,2)./(sum(DOS_M,2))...
         ) *1e10;
            % D_ivs_e_M(i_b,i_b2,1) = ...
            %     sqrt(sum(sum(DP_temp.^2.*DOS_M))/(sum(sum(DOS_M)))) *1e10;

            hbar_w_ivs_e_M(i_b,i_b2,:) = w_v2_array;
        end
    end

end


for i_b = 1:n_VB
    
    for im = n_m:-1:1    
        for iq = n_q:-1:1

        DP_temp(im,iq) = DP_tw(im,iq).VB(i_b,i_b) ;
        
        end
    end

     D_odp_h_M(:,i_b) = ...
         sqrt(...
         sum(DP_temp.^2.*DOS_M,2)./(sum(DOS_M,2))...
         ) *1e10;
         % D_odp_h_M(1,i_b) = ...
         %     sqrt(sum(sum(DP_temp.^2.*DOS_M))/(sum(sum(DOS_M)))) *1e10;

     hbar_w_odp_h_M(:,i_b) = w_v2_array;
            
     for i_b2 = 1:n_VB
        if ne(i_b,i_b2)
            for im = n_m:-1:1    
                for iq = n_q:-1:1
                    DP_temp(im,iq) =  DP_tw(im,iq).VB(i_b,i_b2) ;
                end
            end

            D_ivs_h_M(i_b,i_b2,:) = ...
                sqrt(...
         sum(DP_temp.^2.*DOS_M,2)./(sum(DOS_M,2))...
         ) *1e10;
            % D_ivs_h_M(i_b,i_b2,1) = ...
            %     sqrt(sum(sum(DP_temp.^2.*DOS_M))/(sum(sum(DOS_M)))) *1e10;

            hbar_w_ivs_h_M(i_b,i_b2,:) = w_v2_array;
        end
    end

end

if min(hbar_w_odp_e_M)<0
    hw_f = min(hbar_w_odp_h_M(hbar_w_odp_h_M>0)) / 2;
    hbar_w_odp_e_M(hbar_w_odp_e_M<0) = hw_f;
    hbar_w_odp_h_M(hbar_w_odp_h_M<0) = hw_f;
    hw_f = min(hbar_w_ivs_h_M(hbar_w_ivs_h_M>0)) / 2;
    hbar_w_ivs_e_M(hbar_w_ivs_e_M<0) = hw_f;
    hbar_w_ivs_h_M(hbar_w_ivs_h_M<0) = hw_f;
end

if strcmp(also_intra,'yes')
    clear DP_temp DOS_M 
    fr = DOS_ph(DOS_ph(:,1)>f_cut & DOS_ph(:,1)<f_cut_intra,1);
    dosr = DOS_ph(DOS_ph(:,1)>f_cut & DOS_ph(:,1)<f_cut_intra,2);
    load(intra_filename)

    for i_b = 1:n_CB

        % for n_i_p = length(n_m_intra_array):-1:1

            for im = n_m_intra_array :-1:1    
                for iq = n_q :-1:1
        
                f_t = DP_tw(im,iq).freq;
        
                DOS_t= interp1(fr,dosr,f_t, 'makima');
                if DOS_t < 0
                    DOS_t = 0;
                end
                DOS_M(im,iq) = DOS_t;
        
                DP_temp(im,iq) = DP_tw(im,iq).CB(i_b,i_b) ;
                
                end
            end

            for im = n_m_intra_array:-1:1
                for iq = n_q:-1:1
                    freq_temp_intra(im,iq) = DP_tw(im,iq).freq ;
                end
            end
            w_v2_array_intra = ( (sum(freq_temp_intra.*DOS_M,2))./sum(DOS_M,2) )*4.136e-3 ; % in eV
    
        
             % D_odp_e_M(n_i_p+1,i_b) = ...
             %     sqrt(sum(sum(DP_temp.^2.*DOS_M,2))/(sum(DOS_M,2))) *1e10;
        
             D_odp_e_M(n_m+1:n_m_intra_array+n_m,i_b) = ...
                 sqrt(...
                 sum(DP_temp.^2.*DOS_M,2)./(sum(DOS_M,2))...
                 ) *1e10;
             hbar_w_odp_e_M(n_m+1:n_m_intra_array+n_m,i_b) = w_v2_array_intra;
                
             for i_b2 = 1:n_CB
                if ne(i_b,i_b2)
                    for im = n_m_intra_array:-1:1    
                        for iq = n_q:-1:1
                            DP_temp(im,iq) =  DP_tw(im,iq).CB(i_b,i_b2) ;
                        end
                    end
                    % 
                    % D_ivs_e_M(i_b,i_b2,n_i_p+1) = ...
                    %     sqrt(sum(sum(DP_temp.^2.*DOS_M,2))/(sum(DOS_M,2))) *1e10;
        
                    D_ivs_e_M(i_b,i_b2,n_m+1:n_m+n_m_intra_array) = ...
                            sqrt(...
                     sum(DP_temp.^2.*DOS_M,2)./(sum(DOS_M,2))...
                     ) *1e10;
                    hbar_w_ivs_e_M(i_b,i_b2,n_m+1:n_m+n_m_intra_array) = w_v2_array_intra;
                end
             end

        % end
    
    end
    
    
    for i_b = 1:n_VB

        % for n_i_p = length(n_m_intra_array):-1:1
            
            for im = n_m_intra_array :-1:1 
                for iq = n_q:-1:1
    
                    DP_temp(im,iq) = DP_tw(im,iq).VB(i_b,i_b) ;
                    
                end
            end

            D_odp_h_M(n_m+1:n_m_intra_array+n_m,i_b) = ...
                 sqrt(...
                 sum(DP_temp.^2.*DOS_M,2)./(sum(DOS_M,2))...
                 ) *1e10;
             hbar_w_odp_h_M(n_m+1:n_m_intra_array+n_m,i_b) = w_v2_array_intra;
            
            %  D_odp_h_M(n_i_p+1,i_b) = ...
            %      sqrt(sum(sum(DP_temp.^2.*DOS_M,2))/(sum(DOS_M,2))) *1e10;
            % 
            % hbar_w_odp_h_M(n_i_p+1,i_b) = hbar_w_eff_intra;
                    
             for i_b2 = 1:n_CB
                if ne(i_b,i_b2)
                    for im = n_m_intra_array:-1:1    
                        for iq = n_q:-1:1
                            DP_temp(im,iq) =  DP_tw(im,iq).VB(i_b,i_b2) ;
                        end
                    end

                    D_ivs_h_M(i_b,i_b2,n_m+1:n_m+n_m_intra_array) = ...
                            sqrt(...
                     sum(DP_temp.^2.*DOS_M,2)./(sum(DOS_M,2))...
                     ) *1e10;
                    hbar_w_ivs_h_M(i_b,i_b2,n_m+1:n_m+n_m_intra_array) = w_v2_array_intra;

                    % D_ivs_h_M(i_b,i_b2,n_i_p+1) = ...
                    %     sqrt(sum(sum(DP_temp.^2.*DOS_M,2))/(sum(DOS_M,2))) *1e10;
                    % 
                    % hbar_w_ivs_h_M(i_b,i_b2,n_i_p+1) = hbar_w_eff_intra;
                end
            end
        % end
    end


end


%
k_s = 3;
k_inf = 3; 
Z_i = 1;

save(save_filename, "Z_i", 'k_inf', 'k_s', 'hbar_w_ivs_h_M', 'D_ivs_h_M'...
    , "hbar_w_odp_h_M", "D_odp_h_M", "hbar_w_ivs_e_M", "hbar_w_odp_e_M" ...
    , "D_odp_e_M", "D_ivs_e_M", "rho_mass_density", "bands_type_VB", "bands_type_CB")
