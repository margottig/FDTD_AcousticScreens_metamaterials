clc; clear all; close all;

metalens=["meta1.bmp", "meta2.bmp", "meta3.bmp"];
fcent=[125 250 500 1000 2000 4000];

h3_matrix = [];
href = []; %declaring href for later use
for ii=1:length(metalens)
    tic
    [h,h2,h3]=fdtd2d_metalens(metalens(ii), 2000);
    toc
    href=h;
    h3_matrix = [h3_matrix; h3];
end
metalens = cellstr(metalens);
plot_IL(metalens, h3_matrix, href);