data_file={'Human_20 yo RIEB15-1632_OD_data.xls', 'Human_21 yo RIEB14-1748_OD_data.xls', 'Human_24 yo RIEB13-0253_OS_data', 'Human_25 yo RIEB15-1976_OD_data.xls', 'Human_26 yo RIEB14-1243_OD_data.xls','Human_27 yo RIEB16-0368_OD_data.xls', 'Human_28 yo RIEB13-1936_OD_data.xls', 'Human_29 yo RIEB13-0768T1_OD_data.xls', 'Human_30 yo RIEB13-0161_OD_data.xls'};
zone_vec=1:6;
age=[20:21,24:30];
label_cell=cell(11,11);
label_cell(1,:)={'Anterior Data', 'Chien','Forbes','Fourier', 'ellipse','', 'Posterior Data', 'Chien','Forbes','Fourier', 'ellipse'};
label_cell(2:end,1)={'Fit', 'Bending Energy', 'Waviness', 'Variance of Curvature', 'ARoC 1', 'ARoC 2', 'ARoC 3', 'ARoC 4', 'ARoC 5', 'ARoC 6'};
label_cell(2:end,7)={'Fit', 'Bending Energy', 'Waviness', 'Variance of Curvature', 'PRoC 1', 'PRoC 2', 'PRoC 3', 'PRoC 4', 'PRoC 5', 'PRoC 6'};
for dat_int=1:length(data_file)
    full_cell=label_cell;
data_matrix=anterior(data_file{dat_int},age(dat_int),zone_vec);
full_cell(2:end, 2:5)=num2cell(data_matrix);
data_matrix=posterior(data_file{dat_int},age(dat_int),zone_vec);
full_cell(2:end, 8:11)=num2cell(data_matrix);
    writecell(full_cell,'Curve data.xls','sheet', strcat('age',num2str(age(dat_int))))
end


