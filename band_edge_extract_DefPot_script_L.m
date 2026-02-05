% n_q=6; % 4;
% n_modes=6;
% n_atoms=2;
modes_info
fileName = 'Si_fs.bxsf';
material_name = 'Si';
alat = 0.54;
reduce_flag = 'n'; % 'y' if SOC
interpolation_factor = 0;

% d0 = 0.01e-10; % displacement is defined as 0.01 Angstrom, results in eV/m

shifting_flag = 'no';
[Ek_i_unshifted, ~ ] = bxsf_to_ELECTRA_shifting_option(fileName,material_name,alat,reduce_flag,interpolation_factor,shifting_flag) ;
shifting_flag = 'yes';
[Ek_i_shifted, Fermi]  = bxsf_to_ELECTRA_shifting_option(fileName,material_name,alat,reduce_flag,interpolation_factor,shifting_flag) ;

% definition of the eventual ref. value
Ref_i = max(max(max(Ek_i_unshifted(:,:,:,1)))); % mean(mean(mean(Ek_i_unshifted(:,:,:,1))));
E_threshold = 0.45; % eV

% definition of VB and CB indexes
[CB_idx_temp, VB_idx_temp, ~ ] = shifting_bands(Ek_i_shifted) ; % the Ek with CB edge set to zero does not interest now

for id_band = (max(size(CB_idx_temp))):-1:1
    E_temp(id_band) = min(min(min(Ek_i_shifted(:,:,:,CB_idx_temp(id_band)))));
end
CB_edge = min(E_temp); 
clear E_temp;
CB_idx = 0*CB_idx_temp;
for id_band = (max(size(CB_idx_temp))):-1:1
    if min(min(min(Ek_i_shifted(:,:,:,CB_idx_temp(id_band))))) - CB_edge < 0.5
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
    if abs(max(max(max(Ek_i_shifted(:,:,:,VB_idx_temp(id_band))))) - VB_edge) < 0.5
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
    CB_i(id_band).id_km = [];
    CB_i(id_band).Em = [];
    CB_i(id_band).min_array = [];
    CB_i(id_band).id_k_range = [];
    CB_i(id_band).E_ave = 0;
    CB_i(id_band).E_ave_unshifted = 0;
    CB_i(id_band).E_ave_cut = 0;
    CB_i(id_band).E_ave_unshifted_cut = 0; 
end
CB_Gamma_flag = zeros(1,max(size(CB_idx)));
% evaluation
for id_band = max(size(CB_idx)):-1:1

    Ek_temp = Ek_i_shifted_CB(:,:,:,CB_idx(id_band)) ;
    Ek_temp_unshifted = Ek_i_unshifted(:,:,:,CB_idx(id_band)) ;
    
    id_k = imregionalmin(Ek_temp,26) ;
    Em = Ek_temp(id_k) ;
    min_array = Em(Em < min(min(min(Ek_temp)))+0.05);
    
    CB_i(id_band).id_km = id_k;
    CB_i(id_band).Em = Em;
    CB_i(id_band).min_array = min_array;


    id_k = find( (Ek_temp) < E_threshold );
    E_range = Ek_temp(id_k);
    E_range_unshifted = Ek_temp_unshifted(id_k);

    CB_i(id_band).id_k_range = id_k;
    CB_i(id_band).E_ave = mean(E_range);
    CB_i(id_band).E_ave_unshifted = mean(E_range_unshifted);

end

if strcmp(semiconductor,'yes')
    %initialization
    VB_i = struct();
    for id_band = (max(size(VB_idx))):-1:1
        VB_i(id_band).id_km = [];
        VB_i(id_band).min_array = [];
        VB_i(id_band).id_k_range = [];
        VB_i(id_band).E_ave = 0;
        VB_i(id_band).E_ave_unshifted = 0;
        VB_i(id_band).E_ave_cut = 0;
        VB_i(id_band).E_ave_unshifted_cut = 0;
    end
    VB_Gamma_flag = zeros(1,max(size(CB_idx)));
    % evaluation
    for id_band = max(size(VB_idx)):-1:1
    
        Ek_temp = -Ek_i_shifted_VB(:,:,:,VB_idx(id_band)) ;
        Ek_temp_unshifted = -Ek_i_unshifted(:,:,:,VB_idx(id_band)) ;
        
        id_k = imregionalmin(Ek_temp,26) ;
        Em = Ek_temp(id_k) ;
        min_array = Em(Em < min(min(min(Ek_temp)))+0.05);
    
        VB_i(id_band).id_km = id_k;
        VB_i(id_band).Em = -Em;
        VB_i(id_band).min_array = -min_array;
    
        id_k = find( (Ek_temp) < E_threshold );
        E_range = Ek_temp(id_k);
        E_range_unshifted = Ek_temp_unshifted(id_k);
    
        VB_i(id_band).id_k_range = id_k;
        VB_i(id_band).E_ave = -mean(E_range);
        VB_i(id_band).E_ave_unshifted = -mean(E_range_unshifted);
    
    end
end

% initialization of the DP structures
DP_shifted = struct();
DP_unshifted = struct();
coord_struct = load('coord_modes.mat','coord');
coord = coord_struct.coord;
for im = 1:n_modes
    for iq = 1:n_q
        DP_shifted(im,iq).CB = zeros(max(size(CB_idx)), max(size(CB_idx)));
        DP_shifted(im,iq).CB_ivs = zeros(max(size(CB_idx)), max(size(CB_idx)));
        DP_unshifted(im,iq).CB = zeros(max(size(CB_idx)), max(size(CB_idx)));
        DP_unshifted(im,iq).CB_ivs = zeros(max(size(CB_idx)), max(size(CB_idx)));
        DP_shifted(im,iq).CB_Z = zeros(max(size(CB_idx)), max(size(CB_idx)));
        DP_unshifted(im,iq).CB_Z = zeros(max(size(CB_idx)), max(size(CB_idx)));
        DP_shifted(im,iq).freq = 0;
        DP_unshifted(im,iq).freq = 0;
        if strcmp(semiconductor,'yes')
            DP_shifted(im,iq).VB = zeros(max(size(VB_idx)), max(size(VB_idx)));
            DP_shifted(im,iq).VB_ivs = zeros(max(size(VB_idx)), max(size(VB_idx)));
            DP_unshifted(im,iq).VB = zeros(max(size(VB_idx)), max(size(VB_idx)));
            DP_unshifted(im,iq).VB_ivs = zeros(max(size(VB_idx)), max(size(VB_idx)));
            DP_shifted(im,iq).VB_Z = zeros(max(size(VB_idx)), max(size(VB_idx)));
            DP_unshifted(im,iq).VB_Z = zeros(max(size(VB_idx)), max(size(VB_idx)));
        end
    end
end

for im = 1:n_modes
    for iq = 1:n_q
        path_to_dir = [num2str(im),'/',num2str(iq),'/'];
        cd(path_to_dir) % the bxsf_to... files must be copied in every dir.
        shifting_flag = 'no';
        [Ek_disp_unshifted, ~ ] = bxsf_to_ELECTRA_shifting_option(fileName,material_name,alat,reduce_flag,interpolation_factor,shifting_flag) ;
        shifting_flag = 'yes';
        [Ek_disp_shifted, Fermi]  = bxsf_to_ELECTRA_shifting_option(fileName,material_name,alat,reduce_flag,interpolation_factor,shifting_flag) ;
        
        Ref_disp = max(max(max(Ek_disp_unshifted(:,:,:,1)))); % mean(mean(mean(Ek_disp_unshifted(:,:,:,1))));

        d0 = coord(1,im,iq).ave_displ*1e-10 ;
        % initialize the structure(s)
        CB_disp = struct();
        for id_band = (max(size(CB_idx))):-1:1
            CB_disp(id_band).Em_displaced = [];
            CB_disp(id_band).min_array_displaced = [];
            CB_disp(id_band).E_ave_displaced = 0;
            CB_disp(id_band).E_ave_unshifted_displaced = 0;
            CB_disp(id_band).E_ave_displaced_cut = 0;
            CB_disp(id_band).E_ave_unshifted_displaced_cut = 0;
            CB_disp(id_band).Z = [];
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
                VB_disp(id_band).Em_displaced = 0;
                VB_disp(id_band).min_array_displaced = 0;
                VB_disp(id_band).E_ave_displaced = 0;
                VB_disp(id_band).E_ave_unshifted_displaced = 0;
                VB_disp(id_band).E_ave_displaced_cut = 0;
                VB_disp(id_band).E_ave_unshifted_displaced_cut = 0;
                VB_disp(id_band).Z = 0;
            end

        else
            Ek_disp_shifted_CB = Ek_disp_shifted;
        end


        % the part below, until line 358, must be repeated for VB
        for id_band = max(size(CB_idx)):-1:1
            Ek_temp = Ek_disp_shifted_CB(:,:,:,CB_idx(id_band)) ;
            Ek_temp_unshifted = Ek_disp_unshifted(:,:,:,CB_idx(id_band)) ;
            
            id_k = CB_i(id_band).id_km ;
            Em_displaced = Ek_temp(id_k) ;
            CB_disp(id_band).Em_displaced = Em_displaced;
            CB_disp(id_band).min_array_displaced = Em_displaced(Em_displaced < min(min(min(Ek_temp)))+0.05);
        
            id_k = CB_i(id_band).id_k_range;
            E_range = Ek_temp(id_k);
            E_range_unshifted = Ek_temp_unshifted(id_k);        
            
            CB_disp(id_band).E_ave_displaced = mean(E_range);
            CB_disp(id_band).E_ave_unshifted_displaced = mean(E_range_unshifted);

            if max(size(CB_disp(id_band).min_array_displaced)) ~= max(size(CB_i(id_band).min_array))
                nk = min([max(size(CB_disp(id_band).min_array_displaced)), max(size(CB_i(id_band).min_array))]);
                temp = mink(CB_disp(id_band).min_array_displaced,nk);
                CB_disp(id_band).min_array_displaced = temp;
                temp = mink(CB_i(id_band).min_array,nk);
                CB_i(id_band).min_array = temp;
            end
            if var(CB_disp(id_band).min_array_displaced) ~= var(CB_i(id_band).min_array)
                ivs_flag(id_band) = 1;
                
                id_k = CB_i(id_band).id_km ; % find coordinates of minima
                Ek_filename = ['Ek_',material_name,'.mat'];
                k_matrixes = load(Ek_filename,'kx_matrix','ky_matrix','kz_matrix'); 
                kx_t = k_matrixes.kx_matrix(id_k);
                ky_t = k_matrixes.ky_matrix(id_k);
                kz_t = k_matrixes.kz_matrix(id_k);
                n_m = find(CB_i(id_band).Em < min(min(min(Ek_i_shifted_CB(:,:,:,CB_idx_temp(id_band)))))+0.05);
                kx_m = kx_t(n_m);
                ky_m = ky_t(n_m);
                kz_m = kz_t(n_m);

                diff_minima = CB_i(id_band).min_array - CB_disp(id_band).min_array_displaced ;
                n_m_most = find(abs(diff_minima) == max(abs(diff_minima)));

                id_k = CB_i(id_band).id_k_range; % find coordinates of the n.n. around the "most minima"
                kx_t = k_matrixes.kx_matrix(id_k);
                ky_t = k_matrixes.ky_matrix(id_k);
                kz_t = k_matrixes.kz_matrix(id_k);
                k_cloud = pointCloud([kx_t,ky_t,kz_t]);
                nn_k = floor(max(size(kx_t))/max(size(n_m)));
                id_p_all = zeros(max(size(n_m_most))*nn_k ,1);
                for i = max(size(n_m_most)):-1:1
                    k_point = [kx_m(n_m_most(i)), ky_m(n_m_most(i)), kz_m(n_m_most(i))];                    
                    [id_p, ~] = findNearestNeighbors(k_cloud, k_point, nn_k ) ;
                    id_p_all( nn_k*(i-1)+1:nn_k*i ) = id_p;
                end
                E_ave_displaced = mean(E_range(id_p_all));
                E_ave_unshifted_displaced = mean(E_range_unshifted(id_p_all));
                Ek_i_temp = Ek_i_shifted_CB(:,:,:,CB_idx(id_band)) ;
                E_i_range = Ek_i_temp(id_k);
                E_ave_i= mean(E_i_range(id_p_all));                
                Ek_i_temp = Ek_i_unshifted(:,:,:,CB_idx(id_band)) ;
                E_i_range = Ek_i_temp(id_k);
                E_ave_i_unshifted = mean(E_i_range(id_p_all));

                CB_disp(id_band).E_ave_displaced_cut = E_ave_displaced;
                CB_disp(id_band).E_ave_unshifted_displaced_cut = E_ave_unshifted_displaced;
                CB_i(id_band).E_ave_cut = E_ave_i;
                CB_i(id_band).E_ave_unshifted_cut = E_ave_i_unshifted;


                CB_disp(id_band).Z = max(size(n_m_most)) / max(size(n_m));
            else
                CB_disp(id_band).Z = 1;
                ivs_flag(id_band) = 0;
            end
            n_m = find(CB_i(id_band).Em < min(min(min(Ek_i_shifted_CB(:,:,:,CB_idx_temp(id_band)))))+0.05);
            if max(size(n_m)) == 1 || ( max(size(n_m)) == 8 && ...
                    min(abs(kx_m)) == 0 && min(abs(ky_m)) == 0 && min(abs(kz_m)) == 0 )
                if max(size(n_m)) == 1
                    disp('one only minimum')
                else
                    disp('minimum at Gamma point')
                    CB_Gamma_flag(id_band) = 1;
                end
                CB_disp(id_band).Z = 1;
            end
        end


        % compute hte deformation potentials DP
        for id_band = max(size(CB_idx)):-1:1

            DP_shifted(im,iq).CB(id_band,id_band) = ...
                abs( CB_disp(id_band).E_ave_displaced - CB_i(id_band).E_ave ) / d0 ; % matrixes nb x nb

            DP_unshifted(im,iq).CB(id_band,id_band) = ...
                abs( CB_disp(id_band).E_ave_unshifted_displaced - Ref_disp - ...
                ( CB_i(id_band).E_ave_unshifted - Ref_i ) ) / d0 ; 

            DP_shifted(im,iq).CB_Z(id_band,id_band) = CB_disp(id_band).Z;
            DP_unshifted(im,iq).CB_Z(id_band,id_band) = CB_disp(id_band).Z;

            for id_band_2 = max(size(CB_idx)):-1:1

                DP_shifted(im,iq).CB(id_band,id_band_2) = ...
                0.5*abs( ( CB_disp(id_band).E_ave_displaced + CB_disp(id_band_2).E_ave_displaced ) ...
                - ( CB_i(id_band).E_ave + CB_i(id_band_2).E_ave ) ) / d0 ;
                DP_shifted(im,iq).CB(id_band_2,id_band) = DP_shifted(im,iq).CB(id_band,id_band_2);

                DP_unshifted(im,iq).CB(id_band,id_band_2) = ...
                0.5*abs( (CB_disp(id_band).E_ave_unshifted_displaced - Ref_disp ...
                + CB_disp(id_band_2).E_ave_unshifted_displaced - Ref_disp ) - ...
                ( (CB_i(id_band).E_ave_unshifted - Ref_i) ...
                + (CB_i(id_band_2).E_ave_unshifted - Ref_i) ) ) / d0 ;
                DP_unshifted(im,iq).CB(id_band_2,id_band) = DP_unshifted(im,iq).CB(id_band,id_band_2);

                DP_shifted(im,iq).CB_Z(id_band,id_band_2) = 0.5*(CB_disp(id_band).Z+CB_disp(id_band_2).Z);
                DP_shifted(im,iq).CB_Z(id_band_2,id_band) = DP_shifted(im,iq).CB_Z(id_band,id_band_2);
                DP_unshifted(im,iq).CB_Z(id_band,id_band_2) = DP_shifted(im,iq).CB_Z(id_band,id_band_2);
                DP_unshifted(im,iq).CB_Z(id_band_2,id_band) = DP_shifted(im,iq).CB_Z(id_band,id_band_2);
            end


%             if ivs_flag(id_band) == 1

                DP_shifted(im,iq).CB_ivs(id_band,id_band) = ...
            abs( CB_disp(id_band).E_ave_displaced_cut - CB_i(id_band).E_ave_cut ) / d0 ;

                DP_unshifted(im,iq).CB_ivs(id_band,id_band) = ...
            abs( CB_disp(id_band).E_ave_unshifted_displaced_cut - Ref_disp - ...
            ( CB_i(id_band).E_ave_unshifted_cut - Ref_i ) ) / d0 ;
                
                for id_band_2 = max(size(CB_idx)):-1:1
                    DP_shifted(im,iq).CB_ivs(id_band,id_band_2) = ...
                    0.5*abs( ( CB_disp(id_band).E_ave_displaced_cut + CB_disp(id_band_2).E_ave_displaced_cut ) ...
                    - ( CB_i(id_band).E_ave_cut + CB_i(id_band_2).E_ave_cut ) ) / d0 ;
                    if DP_shifted(im,iq).CB_ivs(id_band,id_band_2) == 0
                        DP_shifted(im,iq).CB_ivs(id_band,id_band_2) = DP_shifted(im,iq).CB(id_band,id_band_2);
                    end
                    DP_shifted(im,iq).CB_ivs(id_band_2,id_band) = DP_shifted(im,iq).CB_ivs(id_band,id_band_2);
    
                    DP_unshifted(im,iq).CB_ivs(id_band,id_band_2) = ...
                    0.5*abs( (CB_disp(id_band).E_ave_unshifted_displaced_cut - Ref_disp ...
                    + CB_disp(id_band_2).E_ave_unshifted_displaced_cut - Ref_disp ) - ...
                    ( (CB_i(id_band).E_ave_unshifted_cut - Ref_i) ...
                    + (CB_i(id_band_2).E_ave_unshifted_cut - Ref_i) ) ) / d0 ;
                    if DP_unshifted(im,iq).CB_ivs(id_band,id_band_2) == 0
                        DP_unshifted(im,iq).CB_ivs(id_band,id_band_2) = DP_unshifted(im,iq).CB(id_band,id_band_2);
                    end
                    DP_unshifted(im,iq).CB_ivs(id_band_2,id_band) = DP_unshifted(im,iq).CB_ivs(id_band,id_band_2);
                end

%             end
        end


        % as above for VB, from line 208
        if strcmp(semiconductor,'yes')
            ivs_flag_CB = ivs_flag;
            clear ivs_flag
            for id_band = max(size(VB_idx)):-1:1
                Ek_temp = -Ek_disp_shifted_VB(:,:,:,VB_idx(id_band)) ;
                Ek_temp_unshifted = -Ek_disp_unshifted(:,:,:,VB_idx(id_band)) ;
                
                id_k = VB_i(id_band).id_km ;
                Em_displaced = Ek_temp(id_k) ;
                VB_disp(id_band).Em_displaced = -Em_displaced;
                VB_disp(id_band).min_array_displaced = -Em_displaced(Em_displaced < min(min(min(Ek_temp)))+0.05);
            
                id_k = VB_i(id_band).id_k_range;
                E_range = Ek_temp(id_k);
                E_range_unshifted = Ek_temp_unshifted(id_k);        
                
                VB_disp(id_band).E_ave_displaced = -mean(E_range);
                VB_disp(id_band).E_ave_unshifted_displaced = -mean(E_range_unshifted);
    
                if max(size(VB_disp(id_band).min_array_displaced)) ~= max(size(VB_i(id_band).min_array))
                    nk = min([max(size(VB_disp(id_band).min_array_displaced)), max(size(VB_i(id_band).min_array))]);
                    temp = mink(-VB_disp(id_band).min_array_displaced,nk);
                    VB_disp(id_band).min_array_displaced = -temp;
                    temp = mink(-VB_i(id_band).min_array,nk);
                    VB_i(id_band).min_array = -temp;
                end
                if var(VB_disp(id_band).min_array_displaced) ~= var(VB_i(id_band).min_array)
                    ivs_flag(id_band) = 1;
                        
                    id_k = VB_i(id_band).id_km ; % find coordinates of minima
                    Ek_filename = ['Ek_',material_name,'.mat'];
                    k_matrixes = load(Ek_filename,'kx_matrix','ky_matrix','kz_matrix'); 
                    kx_t = k_matrixes.kx_matrix(id_k);
                    ky_t = k_matrixes.ky_matrix(id_k);
                    kz_t = k_matrixes.kz_matrix(id_k);
                    n_m = find(VB_i(id_band).Em > min(min(min(Ek_i_shifted_VB(:,:,:,VB_idx_temp(id_band)))))-0.05);
                    kx_m = kx_t(n_m);
                    ky_m = ky_t(n_m);
                    kz_m = kz_t(n_m);
    
                    diff_minima = VB_i(id_band).min_array - VB_disp(id_band).min_array_displaced;
%                     diff_minima = VB_i(id_band).Em - VB_disp(id_band).Em_displaced ;
                    n_m_most = find(abs(diff_minima) == max(abs(diff_minima)));
%                     kx_m = kx_t(n_m_most);
%                     ky_m = ky_t(n_m_most);
%                     kz_m = kz_t(n_m_most);
    
                    id_k = VB_i(id_band).id_k_range; % find coordinates of the n.n. around the "most minima"
                    kx_t = k_matrixes.kx_matrix(id_k);
                    ky_t = k_matrixes.ky_matrix(id_k);
                    kz_t = k_matrixes.kz_matrix(id_k);
                    k_cloud = pointCloud([kx_t,ky_t,kz_t]);
                    nn_k = floor(max(size(kx_t))/max(size(n_m)));
                    id_p_all = zeros(max(size(n_m_most))*nn_k ,1);
                    for i = max(size(n_m_most)):-1:1
                        k_point = [kx_m(n_m_most(i)), ky_m(n_m_most(i)), kz_m(n_m_most(i))];
%                         k_point = [kx_m(i), ky_m(i), kz_m(i)];
                        [id_p, ~] = findNearestNeighbors(k_cloud, k_point, nn_k ) ;
                        id_p_all( nn_k*(i-1)+1:nn_k*i ) = id_p;
                    end
                    E_ave_displaced = mean(E_range(id_p_all));
                    E_ave_unshifted_displaced = mean(E_range_unshifted(id_p_all));
                    Ek_i_temp = -Ek_i_shifted_VB(:,:,:,VB_idx(id_band)) ;
                    E_i_range = Ek_i_temp(id_k);
                    E_ave_i= mean(E_i_range(id_p_all));                
                    Ek_i_temp = Ek_i_unshifted(:,:,:,VB_idx(id_band)) ;
                    E_i_range = Ek_i_temp(id_k);
                    E_ave_i_unshifted = mean(E_i_range(id_p_all));
    
                    VB_disp(id_band).E_ave_displaced_cut = -E_ave_displaced;
                    VB_disp(id_band).E_ave_unshifted_displaced_cut = -E_ave_unshifted_displaced;
                    VB_i(id_band).E_ave_cut = -E_ave_i;
                    VB_i(id_band).E_ave_unshifted_cut = -E_ave_i_unshifted;
    
    
                    VB_disp(id_band).Z = max(size(n_m_most)) / max(size(n_m));
                else
                    VB_disp(id_band).Z = 1;
                    ivs_flag(id_band) = 0;
                end
                n_m = find(VB_i(id_band).Em > min(min(min(Ek_i_shifted_VB(:,:,:,VB_idx_temp(id_band)))))-0.05);
                if max(size(n_m)) == 1 || ( max(size(n_m)) == 8 && ...
                        min(abs(kx_m)) == 0 && min(abs(ky_m)) == 0 && min(abs(kz_m)) == 0 )
                    if max(size(n_m)) == 1
                        disp('one only minimum')
                    else
                        disp('minimum at Gamma point')
                        VB_Gamma_flag(id_band) = 1;
                    end
                    VB_disp(id_band).Z = 1;
                end
            end
    
    
            % compute the deformation potentials DP
            for id_band = max(size(VB_idx)):-1:1
    
                DP_shifted(im,iq).VB(id_band,id_band) = ...
                    abs( VB_disp(id_band).E_ave_displaced - VB_i(id_band).E_ave ) / d0 ; % matrixes nb x nb
    
                DP_unshifted(im,iq).VB(id_band,id_band) = ...
                    abs( VB_disp(id_band).E_ave_unshifted_displaced - Ref_disp - ...
                    ( VB_i(id_band).E_ave_unshifted - Ref_i ) ) / d0 ; 
    
                DP_shifted(im,iq).VB_Z(id_band,id_band) = VB_disp(id_band).Z;
                DP_unshifted(im,iq).VB_Z(id_band,id_band) = VB_disp(id_band).Z;
    
                for id_band_2 = max(size(VB_idx)):-1:1
    
                    DP_shifted(im,iq).VB(id_band,id_band_2) = ...
                    0.5*abs( ( VB_disp(id_band).E_ave_displaced + VB_disp(id_band_2).E_ave_displaced ) ...
                    - ( VB_i(id_band).E_ave + VB_i(id_band_2).E_ave ) ) / d0 ;
                    DP_shifted(im,iq).VB(id_band_2,id_band) = DP_shifted(im,iq).VB(id_band,id_band_2);
    
                    DP_unshifted(im,iq).VB(id_band,id_band_2) = ...
                    0.5*abs( (VB_disp(id_band).E_ave_unshifted_displaced - Ref_disp ...
                    + VB_disp(id_band_2).E_ave_unshifted_displaced - Ref_disp ) - ...
                    ( (VB_i(id_band).E_ave_unshifted - Ref_i) ...
                    + (VB_i(id_band_2).E_ave_unshifted - Ref_i) ) ) / d0 ;
                    DP_unshifted(im,iq).VB(id_band_2,id_band) = DP_unshifted(im,iq).VB(id_band,id_band_2);
    
                    DP_shifted(im,iq).VB_Z(id_band,id_band_2) = 0.5*(VB_disp(id_band).Z+VB_disp(id_band_2).Z);
                    DP_shifted(im,iq).VB_Z(id_band_2,id_band) = DP_shifted(im,iq).VB_Z(id_band,id_band_2);
                    DP_unshifted(im,iq).VB_Z(id_band,id_band_2) = DP_shifted(im,iq).VB_Z(id_band,id_band_2);
                    DP_unshifted(im,iq).VB_Z(id_band_2,id_band) = DP_shifted(im,iq).VB_Z(id_band,id_band_2);
                end
    
    
%                 if ivs_flag(id_band) == 1
    
                    DP_shifted(im,iq).VB_ivs(id_band,id_band) = ...
                abs( VB_disp(id_band).E_ave_displaced_cut - VB_i(id_band).E_ave_cut ) / d0 ;
    
                    DP_unshifted(im,iq).VB_ivs(id_band,id_band) = ...
                abs( VB_disp(id_band).E_ave_unshifted_displaced_cut - Ref_disp - ...
                ( VB_i(id_band).E_ave_unshifted_cut - Ref_i ) ) / d0 ;
                    
                    for id_band_2 = max(size(VB_idx)):-1:1
                        DP_shifted(im,iq).VB_ivs(id_band,id_band_2) = ...
                        0.5*abs( ( VB_disp(id_band).E_ave_displaced_cut + VB_disp(id_band_2).E_ave_displaced_cut ) ...
                        - ( VB_i(id_band).E_ave_cut + VB_i(id_band_2).E_ave_cut ) ) / d0 ;
                        if DP_shifted(im,iq).VB_ivs(id_band,id_band_2) == 0
                            DP_shifted(im,iq).VB_ivs(id_band,id_band_2) = DP_shifted(im,iq).VB(id_band,id_band_2);
                        end
                        DP_shifted(im,iq).VB_ivs(id_band_2,id_band) = DP_shifted(im,iq).VB_ivs(id_band,id_band_2);
        
                        DP_unshifted(im,iq).VB_ivs(id_band,id_band_2) = ...
                        0.5*abs( (VB_disp(id_band).E_ave_unshifted_displaced_cut - Ref_disp ...
                        + VB_disp(id_band_2).E_ave_unshifted_displaced_cut - Ref_disp ) - ...
                        ( (VB_i(id_band).E_ave_unshifted_cut - Ref_i) ...
                        + (VB_i(id_band_2).E_ave_unshifted_cut - Ref_i) ) ) / d0 ;
                        if DP_unshifted(im,iq).VB_ivs(id_band,id_band_2) == 0
                            DP_unshifted(im,iq).VB_ivs(id_band,id_band_2) = DP_unshifted(im,iq).VB(id_band,id_band_2);
                        end
                        DP_unshifted(im,iq).VB_ivs(id_band_2,id_band) = DP_unshifted(im,iq).VB_ivs(id_band,id_band_2);
                    end
    
%                 end
            end        
        
        end


    save mode_point_specific_data
    cd ../../
    DP_shifted(im,iq).freq = coord(1,im,iq).freq;
    DP_unshifted(im,iq).freq = coord(1,im,iq).freq;
    end
end
save('DeformationPotentials_data', "DP_shifted", "DP_unshifted", "Ref_i", "coord", "VB_Gamma_flag", "CB_Gamma_flag")






% 
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