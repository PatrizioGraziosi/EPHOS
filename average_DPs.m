load DeformationPotentials_data.mat

n_m = size(DP_tw,1);
n_q = size(DP_tw,2);
n_b = size(DP_tw(1,1).VB,1);
for im = n_m:-1:1
	for iq = n_q:-1:1
		for ib = n_b:-1:1
			VB_eff(im,iq,ib) = DP_tw(im,iq).VB(ib,ib);
end
                d0_eff(im,iq) = DP_tw(im,iq).d0;
end
end

DP_VB_eff = squeeze(sum(sum(VB_eff,1),2))';
d0_av = mean(mean(d0_eff));

disp(' ')
disp(DP_VB_eff)
disp(' ')
disp (d0_av)
