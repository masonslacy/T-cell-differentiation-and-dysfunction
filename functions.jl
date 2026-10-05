# functions used in main.jl

function gen_scaffold() 
    # Generates scaffold based on experimental measurements, and output the average density of micro-rods in scaffold

    if !isnan(seed_ρ) # set seed if not using default seed
        Random.seed!(seed_ρ); 
    end

    N = M = Int(floor(sqrt(A_well)*10000/h*sqrt(area_scale))); # side length of square simulation domain (based on dish/well size) assuming each grid cell is the size of one T cell

    wid = (M-1)*h; # total domain width and height
    hei = (N-1)*h;

    xpos_PDFM = xpos_PDF(M); # define the probability distribution functions based off the end points in the domain
    ypos_PDFN = ypos_PDF(N);

    # number of rods in the starting domain (input concentration * total area / MSR mass (gives number of rods in total) scaled down by area and mass scales
    num_rods = Int(round(m_input/MSR_m*area_scale*mass_scale));

    # generate domain
    ρ = zeros(N,M); # initialise scaffold domain
    @showprogress 1 "Generating scaffold..." for rod = 1:num_rods # for each rod, place them randomly in the domain
        if use_MvPDF == false # if not using a multivariate PDF
            x0 = rand(xpos_PDFM); y0 = rand(ypos_PDFN); # random starting indices for rod (fractional)
        else 
            x0, y0 = rand(use_MvPDF); # random starting indices for rod (fractional)
        end
        rod_len = max(rand(Normal(mean_MSR_len, sd_MSR_len))/h, 1); # current micro-rod length sampled from normal distribution in terms of lattice site widths (minimum of 1)
        rod_wid = max(rand(Normal(mean_MSR_diam, sd_MSR_diam))/h, 1); # current micro-rod width sampled from normal distribution in terms of lattice site widths (minimum of 1)
        rod_dir = rand(rot_PDF); # current rod direction in grid domain (starting at rod_ind)

        a = tan(rod_dir); b = -1; c = y0-tan(rod_dir)*x0; # equation of line ax+by+c=0 (tan(rod_dir)*ind[2] - ind[1] + y0-tan(rod_dir)*x0 = 0)

        xind_range = Int(round(x0-rod_len*2)):Int(round(x0+rod_len*2)); # range of indices to check over for rod placement
        yind_range = Int(round(y0-rod_len*2)):Int(round(y0+rod_len*2));

        # for all potential cells the rod may fall on
        for x_ind = xind_range
            for y_ind = yind_range 
                dist_n = abs(a*x_ind+b*y_ind+c)/sqrt(a^2+b^2); # distance from new indices to the closest point on the line
                if dist_n <= rod_wid/2 && sqrt((x_ind-x0)^2+(y_ind-y0)^2) <= rod_len/2 # if the point is within half the width and length of the rod
                    x_ind_n = x_ind;
                    if x_ind_n < 1 # 'wrap around' the rods if they exit the domain
                        x_ind_n += M;
                    elseif x_ind_n > M
                        x_ind_n -= M;
                    end
                    y_ind_n = y_ind;
                    if y_ind_n < 1
                        y_ind_n += N;
                    elseif y_ind_n > N
                        y_ind_n -= N;
                    end
                    ρ[y_ind_n,x_ind_n] = min(ρ[y_ind_n,x_ind_n]+1/(rod_len*rod_wid), max_d); # add rod density to this point, ensuring cumulative density does not exceed max_d
                end
            end
        end
    end

    Random.seed!(); # reset seed back to default

    ρ *= MSR_m; # re-formulate ρ as a mass density of micro-rods (μg/h^2)
    ρ *= act_stim; # scale ρ to represent the mass density of activating stimulus on micro-rods (ng/h^2)

    return mean(ρ); # average stimuli mass density
end

function get_data(plot_data)
    # compiles the data from (Zhang, et al., 2023) and computes cell counts based on the data (no data available for IL-2 in this dataset), plots data if plot_data is true
    
    # original data:
    init_T = 1e5; # total number of T cells in the initial sample 


    # initial CD4-to-CD8 ratio (Fig. 1c) [healthy/patient, 10 samples]
    init_CD4toCD8 = [1.305993 1.809386 1.680233 2.614458 5.127716 4.463687 5.133333 5.449275 NaN NaN; 1.28 0.44 0.84 3.87 1.599415 1.67147 2.056738 1.550802 0.34072 0.614828];
    
    # CD4-to-CD8 ratio after expanding 8 days (Fig. S3)
    final_CD4toCD8 = zeros(7,2,10); # [7 stimulus amounts, healthy/patient, 10 samples]
    final_CD4toCD8[:,1,:] = [2.590909091 1.202702703 0.90562249 1.45994832 6.323529412 3.136929461 4.369565217 3.639810427 NaN NaN; 2.521428571 1.178414097 0.80669145 1.050632911 4.655172414 2.662962963 4.52247191 3.180995475 NaN NaN; 2.156549521 1.165577342 0.836158192 2.239202658 3.488038278 2.427046263 4.794117647 4.015544041 NaN NaN; 2.360544218 1.518987342 1.071881607 3.143459916 6.035460993 3.878787879 5.329032258 5.594594595 NaN NaN; 2.85546875 1.943786982 1.451371571 4.37704918 7.684210526 5.64 5.467105263 6.01459854 NaN NaN; 3.245689655 2.155555556 1.693150685 7.31092437 8.631067961 5.775510204 5.043209877 5.957142857 NaN NaN; 4.539325843 2.905511811 1.862973761 6.373134328 9.461457233 6.595419847 5.38961039 6.028776978 NaN NaN];
    final_CD4toCD8[:,2,:] = [1.819484241 0.644444444 0.915354331 NaN NaN 2.61423221 5.736111111 0.774451098 NaN 1.4004914; NaN 1.016736402 1.475949367 NaN NaN 2.343642612 8.08490566 0.480769231 NaN 0.882917466; 3.506849315 1.117256637 1.814492754 6.595419847 NaN 2.424028269 11.69909209 0.373529412 NaN 0.871455577; 2.908730159 1.711048159 1.978593272 NaN NaN 2.670454545 19.85200846 0.619377163 NaN 3.433035714; 3.632075472 2.182724252 2.395104895 NaN NaN 2.76953125 22.35849057 0.670212766 NaN 3.848780488; 4.730994152 2.309027778 2.244147157 NaN NaN 2.942857143 30.44444444 0.921649485 NaN 5.155279503; 4.444444444 2.673076923 2.670454545 NaN NaN 3.511627907 12.58596974 0.858546169 NaN 2.198051948];


    # initial percentage of memory/effector cells (Fig. 1d) 
    init_perc = zeros(4,2,2,10); # [CD45RA and CCR7: ++/-+/--/+-, CD4/CD8, healthy/patient, 10 samples]
    init_perc[1,1,:,:] = [48.57699 46.26627 46 56.9 59.31 63.8 14.8 39.9 NaN NaN; 62.1 11 7.58 58.2 39.7 9.45 46.7 60.6 2.14 37.7]; init_perc[1,2,:,:] = [23.14078 36.26199 60.7 44.9 24.02593 49.9 16.3 26.9 NaN NaN; 32.2 2.01 3.86 35.7 23.6 12.7 15.3 70.4 0.33 33.3]; 
    init_perc[2,1,:,:] = [17.17816 15.48524 15.6 9.46 16.65278 10.9 33.2 20.4 NaN NaN; 6.67 17.5 11.1 15.4 4.98 14.4 25.2 6.1 5.03 15.8]; init_perc[2,2,:,:] = [1.917299 2.221996 2.01 0.55 4.808035 6.19 14.6 10.2 NaN NaN; 0.67 0.99 0.75 4.27 3.01 1.25 2.39 0.11 0.097 1.3]; 
    init_perc[3,1,:,:] = [26.72896 29.21694 32.7 23.6 16.55832 19.6 51.6 37 NaN NaN; 17.5 59.9 58.8 17.9 44.9 53.6 21.1 17.6 77.5 36.7]; init_perc[3,2,:,:] = [23.42884 18.1412 20.7 11.7 38.29333 17.4 44.1 26 NaN NaN; 8.5 36 12.2 7.93 31.6 10.7 11.6 4.22 32.7 18.8]; 
    init_perc[4,1,:,:] = [7.51589 9.031552 5.6 10 7.477635 5.8 2.07 3.93 NaN NaN; 13.8 11.6 22.6 8.51 10.4 22.6 7.05 15.7 15.4 9.86]; init_perc[4,2,:,:] = [51.51309 43.37482 16.6 42.9 32.87271 26.5 25.1 36.9 NaN NaN; 58.7 61 83.2 52.1 41.8 75.4 70.7 25.3 66.9 46.7]; 

    # percentage of memory/effector cells after 8 days (Fig. 1i)
    final_perc = zeros(7,4,2,2,10); # [7 stimulus amounts, CD45RA and CCR7: ++/-+/--/+-, CD4/CD8, healthy/patient, 10 samples]
    final_perc[:,1,1,1,:] = [35.4 33.1 21.6 29.1 22.5 26.5 3.33 7.25 NaN NaN; 27.8 26.5 12.7 14.7 18.7 24.8 5.75 8.5 NaN NaN; 12.3 15.9 13.1 18 22.9 30.9 7.52 10.1 NaN NaN; 6 12.7 14.2 24.7 18.9 36.3 6.22 12.5 NaN NaN; 4.14 12.3 15.3 35.3 15.3 33.7 5.58 13.3 NaN NaN; 3.36 10 16.9 29.9 12.7 30.6 4.93 12.8 NaN NaN; 5.85 11.7 17.3 40.1 12.7 29.5 4.35 10.7 NaN NaN]; final_perc[:,1,1,2,:] = [24.4 6.07 4.85 NaN NaN 7.5 4.72 25.4 NaN 36.3; NaN 6.77 6.78 NaN NaN 8.35 4.16 19.2 NaN 46.4; 14.8 9.41 10 4 NaN 13.6 5.23 22.7 NaN 35; 18.4 13.2 9.69 NaN NaN 19.8 5.48 13.3 NaN 26.1; 16.7 12.4 10.2 NaN NaN 16.9 4.27 12.2 NaN 26.9; 16 11 8.99 NaN NaN 17.5 3.52 12 NaN 31.7; 11.7 8.6 9.33 NaN NaN 13.4 3.43 13 NaN 29.3];
    final_perc[:,2,1,1,:] = [30.6 27.2 26.3 16.7 30.2 15.3 38.1 26.7 NaN NaN; 31.8 23.3 14.2 9.17 21.8 12.1 41 25.6 NaN NaN; 27.9 16.8 10.4 7.52 38.5 17.2 62.5 33.8 NaN NaN; 22.7 16.3 11.4 8.88 48 23.8 67.5 49.4 NaN NaN; 18 15.6 10.6 10.9 54.3 27.8 65.9 49.4 NaN NaN; 17.2 17 12.5 7.09 54.8 34.5 69.8 53.4 NaN NaN; 18.7 14.8 12.9 14.8 54.7 32.5 67.4 54.2 NaN NaN]; final_perc[:,2,1,2,:] = [12.7 24.3 15.6 NaN NaN 17.7 38.2 40.9 NaN 6.9; NaN 28.4 25.9 NaN NaN 21.2 32.9 36.8 NaN 6.19; 26.1 32.3 38.7 23.3 NaN 24.4 36.7 32.3 NaN 5.39; 23.8 45.8 50.6 NaN NaN 46.1 53.8 26.6 NaN 9.67; 33.1 55.1 59.1 NaN NaN 47.5 51.2 25.1 NaN 12.3; 35 59.3 61.1 NaN NaN 56.8 56.7 30.1 NaN 23.3; 34.2 63.4 62.8 NaN NaN 56.3 54.1 28.3 NaN 24.8];
    final_perc[:,3,1,1,:] = [16.3 23.9 36.5 26.8 22.1 24.8 49.2 46.3 NaN NaN; 20.9 23 39.7 24.6 22.9 16.4 41.3 46.2 NaN NaN; 29.3 18.6 25 13.6 20 14 22.6 41.4 NaN NaN; 39.9 21.1 20.3 12.4 20.2 13.1 21 30.1 NaN NaN; 44.9 21 17.3 10.5 20.6 14.2 23.2 29.7 NaN NaN; 50.9 24.2 14.6 8.36 23.6 16.2 20.7 27.5 NaN NaN; 49.4 24.4 12 11.8 23.8 17.2 24.2 29.1 NaN NaN]; final_perc[:,3,1,2,:] = [15.1 48.8 57.1 NaN NaN 36.5 47.6 21.9 NaN 21.1; NaN 46.2 46.5 NaN NaN 31.3 50.7 23.4 NaN 12.9; 23 41.7 35.9 21.7 NaN 21.7 45.9 20.1 NaN 9.22; 20.4 27.6 30.6 NaN NaN 12.9 27.1 17 NaN 6.8; 24.9 22 22.4 NaN NaN 14.5 32.7 18.4 NaN 8.44; 27.3 21.1 22.3 NaN NaN 11.4 31.4 20.8 NaN 9.27; 36.4 21.4 20.1 NaN NaN 15.7 35.2 23.1 NaN 12.1];
    final_perc[:,4,1,1,:] = [17.7 15.7 15.6 27.4 25.2 33.4 9.35 19.7 NaN NaN; 19.5 27.2 33.3 51.5 36.7 46.7 11.9 19.6 NaN NaN; 30.5 48.7 51.5 60.9 18.6 37.9 7.31 14.6 NaN NaN; 31.4 49.9 54.1 54.1 13 26.8 5.28 8.06 NaN NaN; 33 51.1 56.9 43.2 9.82 24.3 5.27 7.52 NaN NaN; 28.6 48.8 56 54.6 8.85 18.7 4.65 6.3 NaN NaN; 26 49.1 57.8 33.3 8.74 20.9 4.02 6 NaN NaN]; final_perc[:,4,1,2,:] = [47.8 20.8 22.5 NaN NaN 38.3 9.44 11.8 NaN 35.7; NaN 18.6 20.8 NaN NaN 39.1 12.3 20.5 NaN 34.5; 36 16.6 15.3 33 NaN 40.3 12.2 24.9 NaN 50.4; 37.4 13.3 9.11 NaN NaN 21.3 13.6 43.1 NaN 57.5; 25.3 10.6 8.26 NaN NaN 21 11.8 44.3 NaN 52.3; 21.8 8.68 7.61 NaN NaN 14.3 8.39 37.1 NaN 35.6; 17.7 6.57 7.81 NaN NaN 14.6 7.33 35.5 NaN 33.8];

    final_perc[:,1,2,1,:] = [41.5 46.6 30.3 26 36 41.9 9.35 8.44 NaN NaN; 30.8 28.9 13 21.7 31 32 14.6 13.2 NaN NaN; 11.5 15.7 16.4 52.6 48.7 51.4 16.2 15.5 NaN NaN; 8.66 23.8 26.7 66.1 42.1 47.7 12.5 13.5 NaN NaN; 12 29.5 33.1 65.6 39.1 44.6 11.4 13.9 NaN NaN; 11.3 23.4 33.1 64 31 39.7 10.6 12.3 NaN NaN; 17.7 26.5 28.1 41.7 33.1 34.1 8.57 9.98 NaN NaN]; final_perc[:,1,2,2,:] = [22.7 6.05 6.37 NaN NaN 18 7.44 34.4 NaN 23.4; NaN 6.59 10.6 NaN NaN 22.6 4.23 19.8 NaN 27.6; 30.6 7.85 16.4 37 NaN 41.9 7.36 19.5 NaN 21.1; 34.9 11.4 16.5 NaN NaN 40.3 12.7 17.3 NaN 50.2; 32.7 11.1 16.5 NaN NaN 35.4 6 14.1 NaN 50; 32.2 9.19 14.4 NaN NaN 30.4 8.79 14.3 NaN 60.2; 26.1 7.58 13.4 NaN NaN 23.9 8.96 13.6 NaN 50.4];
    final_perc[:,2,2,1,:] = [10.1 11.9 16.7 10.1 13.6 11.8 35.9 21.1 NaN NaN; 11.9 7.64 10.7 5.13 7.76 6.93 34.1 21.8 NaN NaN; 5.82 3.5 5.8 5.07 18.3 9.85 49.9 27.1 NaN NaN; 4.58 4.2 5.24 6.13 21.2 12.5 49.6 34.2 NaN NaN; 5.71 4.31 4.46 6.36 24.2 14.7 45.1 29.7 NaN NaN; 7.13 5.08 5.09 5.15 24.8 18.1 42.7 29.9 NaN NaN; 10.5 6.73 4.78 6.08 25.7 16.6 45.5 29.1 NaN NaN]; final_perc[:,2,2,2,:] = [6.99 30.9 19.7 NaN NaN 20.5 24.3 17.9 NaN 2.54; NaN 33.1 32.2 NaN NaN 18.7 22.1 11.2 NaN 1.93; 11 37.9 41.2 11.5 NaN 19.7 20.7 6.14 NaN 2.73; 9.07 46.5 45.4 NaN NaN 29.7 45.2 5.76 NaN 7.18; 15.8 53.4 49.2 NaN NaN 29 39.7 4.76 NaN 7.49; 14.4 52.4 48.2 NaN NaN 36.2 40.8 7.3 NaN 12.1; 14.7 57.9 47.3 NaN NaN 34.1 38 7.38 NaN 11.2];
    final_perc[:,3,2,1,:] = [9.79 11.2 18.5 11.2 12.4 15.4 36.9 43.1 NaN NaN; 11.9 6.39 16.6 7.12 8.31 8 29.9 37.1 NaN NaN; 15 4.66 9.68 4.32 9.47 7.79 23.7 40.4 NaN NaN; 19.1 5.54 8.59 4.48 14.1 10.7 28.4 40.8 NaN NaN; 19.1 5.79 7.74 4.96 16.4 12.4 32.6 43.4 NaN NaN; 25.3 8.66 7.56 3.63 22.4 16 35.3 44.3 NaN NaN; 23.2 9.78 7.78 10.2 21.2 18 36.2 48.7 NaN NaN]; final_perc[:,3,2,2,:] = [5.84 42.5 45.4 NaN NaN 12.8 43.3 14 NaN 8.38; NaN 42.8 35.2 NaN NaN 11.5 44.2 14.2 NaN 5.43; 13.7 41.3 27.8 10.1 NaN 6.67 38.6 8.42 NaN 5.58; 11.9 31 26.7 NaN NaN 8.96 24.3 8.9 NaN 4.08; 17.1 27.1 23.1 NaN NaN 11.6 40.1 10.2 NaN 5.52; 17.9 30.8 26.3 NaN NaN 14.6 32.7 13.9 NaN 5.76; 27.2 28.9 27.8 NaN NaN 20 37.2 18 NaN 7.71];
    final_perc[:,4,2,1,:] = [38.6 30.3 34.5 52.8 38 30.9 17.8 27.4 NaN NaN; 45.3 57.1 59.7 66 52.9 53.1 21.3 27.9 NaN NaN; 67.6 76.1 68.1 38.1 23.6 30.9 10.2 16.9 NaN NaN; 67.7 66.4 59.5 23.3 22.5 29 9.6 11.5 NaN NaN; 63.2 60.4 54.7 23.1 20.3 28.2 10.9 13 NaN NaN; 56.2 62.9 54.3 27.2 21.9 26.2 11.4 13.5 NaN NaN; 48.7 57 59.3 42 20 31.3 9.67 12.2 NaN NaN]; final_perc[:,4,2,2,:] = [64.5 20.5 28.5 NaN NaN 48.7 24.9 33.7 NaN 65.7; NaN 17.4 22 NaN NaN 47.2 29.5 54.8 NaN 65; 44.6 13 14.6 41.4 NaN 31.6 33.3 66 NaN 70.6; 44.2 11.1 11.4 NaN NaN 21.1 17.8 68.1 NaN 38.5; 34.4 8.31 11.1 NaN NaN 24.1 14.2 70.9 NaN 37; 35.4 7.61 11.1 NaN NaN 18.7 17.7 64.5 NaN 21.9; 32 5.63 11.5 NaN NaN 22 15.8 60.9 NaN 30.7];
    

    # initial percentage of CD45RA+CCR7+ T cells that are also CD95+ (not naive) (Fig. S2f) (data digitised with WebPlotDigitizer, roughly ordered the same as other experimental data, but I'm not sure it's the exact same)
    init_CD95 = zeros(2,2,10); # [CD4/CD8, healthy/patient, 10 samples]
    init_CD95[:,1,:] = [13.84440317 21.36190057 5.965176094 4.308727519 8.258430284 3.072802348 33.59818286 7.086234376 NaN NaN; 53.18164576 57.41201241 14.02214341 13.75637972 40.6834976 6.951927962 66.90903473 17.64157666 NaN NaN];
    init_CD95[:,2,:] = [19.41282067 33.59818286 85.80606886 8.205922742 56.32387432 39.91241744 21.0911239 8.153749046 71.77210415 11.07442866; 17.30721333 35.81101418 71.77210415 6.439679688 43.91969139 43.91969139 35.1322829 1.282428736 62.37547278 6.820167076];


    # initial percentage of CD25+ cells (activated)
    init_CD25 = zeros(2,2,10); # [CD4/CD8, healthy/patient, 10 samples]
    init_CD25[:,1,:] = [3.63 2.41 1.01 0.82 3.03 2.95 9.66 2.98 NaN NaN; 0.046 0.065 0.034 0.029 0.17 0.17 0.97 0.11 NaN NaN];
    init_CD25[:,2,:] = [1.53 8.58 9.7 0.63 1.78 2.87 0.43 4.13 9.17 1.43; 0.069 0.16 0.35 0.3 0.065 0.14 0.13 0.061 0.075 0.054];

    # percentage of CD25+ cells (activated) at day 8
    final_CD25 = zeros(7,2,2,10); # [7 stimulus amounts, CD4/CD8, healthy/patient, 10 samples]
    final_CD25[:,1,1,:] = [86.8 81.4 75.6 73.6 95.5 95.1 83 84.9 NaN NaN; 92.6 82.8 76.6 69.6 86.9 84.7 91.6 92 NaN NaN; 92.4 87.9 87.6 85.2 96.9 95.9 98.6 97.1 NaN NaN; 90.4 85.5 90.2 89.8 98.3 99.1 99.4 99.2 NaN NaN; 91.2 87.1 94.2 95.9 98.9 99.3 99.5 99.5 NaN NaN; 90.2 89 96.8 97.7 99.3 99.7 99.6 99.6 NaN NaN; 94.9 94.7 98.4 99.3 99.3 99.7 99.4 99.6 NaN NaN]; final_CD25[:,1,2,:] = [83.7 84.5 76.1 NaN NaN 79.4 95 89.3 NaN 65; NaN 90.2 88.5 NaN NaN 83.4 96.5 91.6 NaN 73.9; 95.1 95.8 95.9 96.6 NaN 89.6 97.7 92.7 NaN 86.5; 96.1 98.9 97.6 NaN NaN 97.6 99.4 97.9 NaN 97.1; 97.9 99.4 99.2 NaN NaN 98.4 99.3 98 NaN 98.2; 99 99.8 99.4 NaN NaN 99.1 99.5 98.8 NaN 99.4; 99.4 99.8 99.7 NaN NaN 99.3 99.5 99 NaN 99.4];
    final_CD25[:,2,1,:] = [67 48.7 56.5 64.8 88 95.3 83.7 79.5 NaN NaN; 68.7 53.6 71.8 79.8 81.9 89.5 92.3 91.8 NaN NaN; 75 81.6 89.3 95.7 96.9 98.3 97.8 96.3 NaN NaN; 78.3 88.5 94 98.4 97.9 99.4 98 98.1 NaN NaN; 83.9 92 96.7 99.3 98.7 99.5 98.2 98.6 NaN NaN; 84.8 93.3 97.6 99.5 99 99.7 98.9 98.8 NaN NaN; 92.3 96.3 98.4 99.5 99.1 99.6 99.1 98.8 NaN NaN]; final_CD25[:,2,2,:] = [66.7 92.6 81.1 NaN NaN 82.3 89.8 70.9 NaN 46.4; NaN 95.9 93.4 NaN NaN 89.6 90.7 75 NaN 49.7; 91.1 98.7 98.5 94.3 NaN 96.3 92.8 77.9 NaN 73.8; 94.2 99.8 98.7 NaN NaN 98.5 99.5 92 NaN 96.5; 95.1 99.8 99.4 NaN NaN 98.7 99.3 91.4 NaN 98; 98.1 99.9 99.6 NaN NaN 98.9 99.5 95.7 NaN 99.4; 99.3 99.8 99.7 NaN NaN 99.3 99.7 97 NaN 99.5];


    # percentage progentior exhausted (PD-1+TIM-3-) at day 0 
    init_TPEX = zeros(2,2,10); # [CD4/CD8, healthy/patient, 10 samples]
    init_TPEX[:,1,:] = [24.5 26.5 20.9 16.1 16.1 9.66 38.2 20.7 NaN NaN; 39.6 28 17.9 24 56.7 30.1 37.5 34.5 NaN NaN];
    init_TPEX[:,2,:] = [12.3 26.7 64.8 10 35.1 30.8 7.02 19.7 50.9 16.6; 38.6 32.5 20.9 13.4 24 25.7 18.5 10.4 20.2 11.6];

    # percentage progentior exhausted (PD-1+TIM-3-) at day 8
    final_TPEX = zeros(7,2,2,10); # [7 stimulus amounts, CD4/CD8, healthy/patient, 10 samples]
    final_TPEX[:,1,1,:] = [15.1 14.9 22.2 11.7 10.5 5.73 35.1 19.6 NaN NaN; 11.6 13.1 17.9 5.96 8.83 5.84 27.7 16.7 NaN NaN; 17.4 15.1 16.3 5.82 9.93 5.18 32.3 16.3 NaN NaN; 19.8 16.1 17.6 7.38 12 5.95 32.3 16.9 NaN NaN; 18 18.6 18.6 10.4 15.4 6.43 31.3 16.3 NaN NaN; 19.1 19.9 16.2 9.05 15.7 7.72 29.2 14.4 NaN NaN; 13.6 17.6 14 17.8 14.4 7.76 32.3 17.2 NaN NaN]; final_TPEX[:,1,2,:] = [7.56 10.1 20.4 NaN NaN 13.1 24.1 26.1 NaN 5.46; NaN 10.8 17.7 NaN NaN 13.6 23.3 21.9 NaN 4.41; 7.45 9.31 20.2 4.56 NaN 12.7 22.6 22.2 NaN 2.45; 5.49 13 24.6 NaN NaN 22.9 22.2 21.2 NaN 1.55; 5.85 17.6 29.9 NaN NaN 23.1 21.2 18.9 NaN 1.89; 8.76 23 32.1 NaN NaN 25.6 22.4 21.8 NaN 4.21; 8.04 23.4 30.3 NaN NaN 27.6 20.6 15.2 NaN 3.27];
    final_TPEX[:,2,1,:] = [3.34 3.68 4.76 4.29 2.65 1.57 9.97 9.35 NaN NaN; 1.58 1.84 2.65 1.82 2.23 1.93 6.41 6.13 NaN NaN; 0.75 1.08 1.53 2.87 2.61 1.97 7.94 5.77 NaN NaN; 1.03 1.52 1.74 4.38 3.16 2.18 6.7 5.26 NaN NaN; 1.55 1.99 2.14 5.32 4.16 2.43 7.07 5.23 NaN NaN; 3.16 2.1 1.83 5.2 4.08 2.68 6.33 4.13 NaN NaN; 3.49 2.87 1.76 8.18 3.71 3.23 7.31 5.75 NaN NaN]; final_TPEX[:,2,2,:] = [0.7 1.78 1.86 NaN NaN 3.01 4.97 2.62 NaN 0.43; NaN 1.65 2.55 NaN NaN 2.11 5.93 1.09 NaN 0.2; 0.55 1.74 3.28 0.68 NaN 2.62 5.63 0.72 NaN 0.11; 0.68 2.21 3.76 NaN NaN 4.12 6.5 0.71 NaN 0.17; 0.8 2.92 5.07 NaN NaN 4.37 12.1 0.58 NaN 0.21; 1.24 3.5 4.96 NaN NaN 4.42 8.29 0.95 NaN 0.52; 1.18 4 5.05 NaN NaN 6.15 7.33 0.74 NaN 0.45];


    # percentage terminally exhausted (PD-1+TIM-3+) at day 0 
    init_TEX = zeros(2,2,10); # [CD4/CD8, healthy/patient, 10 samples]
    init_TEX[:,1,:] = [0.89 0.39 0.12 0.17 0.38 0.33 0.75 0.51 NaN NaN; 0.68 0.52 0.077 0.13 0.8 0.25 0.56 0.45 NaN NaN];
    init_TEX[:,2,:] = [0.35 0.6 1.31 0.63 0.082 0.15 0.033 0.16 0.53 0.2; 0.99 0.51 1.08 2.44 0.16 0.14 0 0.051 0.14 0.096];

    # percentage terminally exhausted (PD-1+TIM-3+) at day 8 (Fig. 1h)
    final_TEX = zeros(7,2,2,10); # [7 stimulus amounts, CD4/CD8, healthy/patient, 10 samples]
    final_TEX[:,1,1,:] = [10.1 10.6 14.9 14.1 32.4 29.2 25 36 NaN NaN; 11.7 13.5 18.8 11.6 23.2 23.3 26 33.9 NaN NaN; 28.1 25.6 27.1 14.2 44.8 40.9 30.6 42.8 NaN NaN; 34.3 28.4 29.2 17.5 53.3 57.4 31.7 54.5 NaN NaN; 40.6 38 37.8 34.2 52.3 59.5 34.3 59.3 NaN NaN; 44.8 47.3 50.2 40 58.3 67.1 40.3 64.3 NaN NaN; 47.9 55.2 60.8 60.7 60.7 69.3 36.7 61.5 NaN NaN]; final_TEX[:,1,2,:] = [17.5 15.7 19 NaN NaN 10.7 16.9 10.5 NaN 19.4; NaN 16.9 22.8 NaN NaN 13.8 24.4 16.8 NaN 20.4; 51.2 24 30.2 34.3 NaN 17.8 32.7 27.6 NaN 20; 65.8 33.2 34.2 NaN NaN 21.4 42.1 55.3 NaN 17.2; 74 34.1 35 NaN NaN 22.7 50.2 58.6 NaN 27.1; 72.6 34 37 NaN NaN 23.6 51.2 57.6 NaN 41.2; 81.4 36.2 41.5 NaN NaN 26.6 58.8 71.2 NaN 49.4];
    final_TEX[:,2,1,:] = [7.12 6.36 6.54 10.1 15.4 23.7 13.4 26.9 NaN NaN; 5.91 4.9 6.71 6.12 7.64 17.6 12.3 23 NaN NaN; 6.95 10 11.5 15.1 21.4 37.3 19.8 32 NaN NaN; 8.94 15.5 19.1 27.9 27 50.1 21.8 40.3 NaN NaN; 12.4 23.1 28 45.8 26.3 46.6 21.9 43.2 NaN NaN; 15.6 28 34.7 48.1 30.9 54.5 30.1 44.7 NaN NaN; 15.9 33.6 37.2 64.8 32.5 59 29.3 45.6 NaN NaN]; final_TEX[:,2,2,:] = [3.29 5.63 4.94 NaN NaN 4.87 14.4 6.98 NaN 9.54; NaN 7.07 7.18 NaN NaN 5.54 24.4 8.01 NaN 5.89; 15.2 9.81 14.6 16.1 NaN 8.92 31.2 10.7 NaN 3.28; 30.3 15.4 21.4 NaN NaN 19.3 39.1 18.3 NaN 5.91; 37.9 17.5 24.1 NaN NaN 18.7 50 17.6 NaN 9.68; 44.6 19.3 27 NaN NaN 21.4 43.4 21.7 NaN 18.7; 60.1 21 31.9 NaN NaN 24.4 58.2 29.1 NaN 22.7];


    # fold expansion at day 8 (Fig. 1g) 
    fold_exp = zeros(7,2,10); # [7 stimulus amounts, healthy/patient, 10 samples]
    fold_exp[:,1,:] = [1.438666667 3.6875 2.7258 6.608 176 215 56.33333333 147.6666667 NaN NaN; 18.48 45.961 99.12 167.56 373.3333333 343.3333333 77 171 NaN NaN; 123.84 267.86 244.85 325.68 416.6666667 396.6666667 105.3333333 277.6666667 NaN NaN; 169.0555556 273.76 290.87 339.25 410 406 138 280.8466667 NaN NaN; 212.6666667 294.41 295 345.74 406.6666667 401.2 179.3333333 302.6666667 NaN NaN; 263.7333333 279.07 241.31 318.6 396.6666667 401.8 104.3333333 300.3333333 NaN NaN; 298.8444444 286.74 254.29 331.58 386.6666667 343.3333333 85.66666667 292.3333333 NaN NaN];
    fold_exp[:,2,:] = [31.4 44 91.5 NaN NaN 116 3.09 4.98 NaN 84; NaN 126 261.5 NaN NaN 200.5 5.82 11.625 NaN 100.5; 41.25 148.5 277.5 200 NaN 260.7 4.41 18.45 NaN 190.5; 57 151.2 231.5 NaN NaN 289.8 11.7 45.45 NaN 175; 47.25 113.22 231.5 NaN NaN 270 9.22 52.95 NaN 169.75; 48.25 112.7 226.5 NaN NaN 273.5 6.71 41.1 NaN 74.75; 30.6 77 208 NaN NaN 164.5 4.59 26.1 NaN 76.5];




    # computations using data: 
    # proportion of day 0 cells that are CD4/CD8 
    init_CD4or8 = zeros(2,2,10); # [CD4/CD8, healthy/patient, 10 samples]
    init_CD4or8[1,:,:] = 1 .- 1 ./ (init_CD4toCD8 .+ 1); # convert ratios to proportions
    init_CD4or8[2,:,:] = 1 ./ (init_CD4toCD8 .+ 1);

    # proportion of day 8 cells that are CD4/CD8 
    final_CD4or8 = zeros(7,2,2,10); # [7 stimulus amounts, CD4/CD8, healthy/patient, 10 samples]
    final_CD4or8[:,1,:,:] = 1 .- 1 ./ (final_CD4toCD8 .+ 1); # convert ratios to proportions
    final_CD4or8[:,2,:,:] = 1 ./ (final_CD4toCD8 .+ 1);

    # average initial percentage of CD45RA+CCR7+ cells that are also CD95+
    mean_CD95 = zeros(2,2); # [CD4/CD8, healthy/patient]
    mean_CD95[:,1] = mean(init_CD95[:,1,1:8],dims=2); # ignore the NaN entries
    mean_CD95[:,2] = mean(init_CD95[:,2,:],dims=2);
    # reformat for computing final number 
    #mean_CD95_comp = zeros(7,2,2);
    #for i = 1:7
    #    mean_CD95_comp[i,:,:] = mean_CD95;
    #end

    # reformat fold expansion data for ease of computation:
    fold_exp_comp = zeros(6,7,2,2,10); # [7 stimulus amounts, CD4/CD8, healthy/patient, 10 samples]
    for cell_type = 1:6
        for CD4or8 = 1:2
            fold_exp_comp[cell_type,:,CD4or8,:,:] = fold_exp;
        end
    end


    
    # ASSUMPTIONS: 
    #  - the amounts of activated and progenitor exhausted are not uniquely identifiable, but we can still compute how many we might expect
    #  - the average percentage of CD45RA+CCR7+ that also express CD95 across all samples is accurate enough to use for individual samples (since we don't know which values correspond to each sample)
    #  - the proportion of cells that are CD45RA+CCR7+CD25- ("non-activated stem cell-like memory") is a good enough proxy for naivety at day 8 (we don't have CD95 expression at day 8)
    #  - activated cells are any cells that express CD25 (consider these separately to naive/memory/effector/TEX - i.e. the sum of those cell types and activated may be over 100% of the cells)
    #  - progenitor exhausted cells express PD-1 but not TIM-3 (consider these separately to naive/memory/effector - i.e. the sum of those cell types and TPEX may be over 100% of the cells)
    #  - TEX cells are any cells (with any expressions) that also express PD-1 and TIM-3 - assuming that these cells are different enough from the other cell types such that a cell may not be memory and terminally exhausted for example


    # initial proportions of each T cell subtype 
    init_prop = zeros(6,2,2,10); # [naive/activated/memory/effector/TPEX/TEX, CD4/CD8, healthy/patient, 10 samples]

    # naive cells: CD45RA+CCR7+CD95-
    init_prop[1,:,:,:] = init_CD4or8 .* init_perc[1,:,:,:]/100 .* (1 .- mean_CD95/100) .* (1 .- init_TEX/100); # the proportion of those that are CD4/CD8 * the proportion of those that are CD45RA+CCR7+ * the proportion of those that are CD95- * the proportion of those that are not terminally exhausted (PD-1+TIM-3+)
    # activated cells: CD25+
    init_prop[2,:,:,:] = init_CD4or8 .* init_CD25/100; # the proportion of those that are CD4/CD8 * the proportion of those that are CD25+
    # memory cells: CD45RA+CCR7+CD95+ and CD45RA-CCR7+
    init_prop[3,:,:,:] = init_CD4or8 .* (init_perc[1,:,:,:]/100 .* mean_CD95/100 + init_perc[2,:,:,:]/100) .* (1 .- init_TEX/100); # the proportion of those that are CD4/CD8 * the proportion of those that are CD45RA+CCR7+ and CD95+ or just CD45RA-CCR7+ * the proportion of those that are not terminally exhausted (PD-1+TIM-3+)
    # effector cells: CD45RA-CCR7- and CD45RA+CCR7-
    init_prop[4,:,:,:] = init_CD4or8 .* (init_perc[3,:,:,:] + init_perc[4,:,:,:])/100 .* (1 .- init_TEX/100); # the proportion of those that are CD4/CD8 * the proportion of those that are CD45RA-CCR7- or CD45RA+CCR7- * the proportion of those that are not terminally exhausted (PD-1+TIM-3+)
    # progenitor exhausted cells: PD-1+TIM-3-
    init_prop[5,:,:,:] = init_CD4or8 .* init_TPEX/100; # the proportion of those that are CD4/CD8 * the proportion of those that are PD-1+TIM-3-
    # terminally exhausted cells: PD-1+TIM-3+
    init_prop[6,:,:,:] = init_CD4or8 .* init_TEX/100; # the proportion of those that are CD4/CD8 * the proportion of those that are PD-1+TIM-3+

    # initial amounts of each T cell subtype
    init_num = init_T * init_prop; # [naive/activated/memory/effector/TPEX/TEX, CD4/CD8, healthy/patient, 10 samples]
    

    # final (day 8) proportions of each T cell subtype 
    final_prop = zeros(6,7,2,2,10); # [naive/activated/memory/effector/TPEX/TEX, 7 stimulus amounts, CD4/CD8, healthy/patient, 10 samples]

    # naive cells: CD45RA+CCR7+CD95-
    final_prop[1,:,:,:,:] = final_CD4or8 .*  final_perc[:,1,:,:,:]/100 .* (1 .- final_CD25/100) .* (1 .- final_TEX/100); # the proportion of those that are CD4/CD8 * the proportion of those that are CD45RA+CCR7+ * the proportion of those that are CD25- * the proportion of those that are not terminally exhausted (PD-1+TIM-3+)
    # activated cells: CD25+
    final_prop[2,:,:,:,:] = final_CD4or8 .* final_CD25/100; # the proportion of those that are CD4/CD8 * the proportion of those that are CD25+
    # memory cells: CD45RA+CCR7+CD95+ and CD45RA-CCR7+
    final_prop[3,:,:,:,:] = final_CD4or8 .* ( final_perc[:,1,:,:,:]/100 .* final_CD25/100 +  final_perc[:,2,:,:,:]/100) .* (1 .- final_TEX/100); # the proportion of those that are CD4/CD8 * the proportion of those that are CD45RA+CCR7+ and CD25+ (not naive) or just CD45RA-CCR7+ * the proportion of those that are not terminally exhausted (PD-1+TIM-3+)
    # effector cells: CD45RA-CCR7- and CD45RA+CCR7-
    final_prop[4,:,:,:,:] = final_CD4or8 .* ( final_perc[:,3,:,:,:] +  final_perc[:,4,:,:,:])/100 .* (1 .- final_TEX/100); # the proportion of those that are CD4/CD8 * the proportion of those that are CD45RA-CCR7- or CD45RA+CCR7- * the proportion of those that are not terminally exhausted (PD-1+TIM-3+)
    # progenitor exhausted cells: PD-1+TIM-3-
    final_prop[5,:,:,:,:] = final_CD4or8 .* final_TPEX/100; # the proportion of those that are CD4/CD8 * the proportion of those that are PD-1+TIM-3-
    # terminally exhausted cells: PD-1+TIM-3+
    final_prop[6,:,:,:,:] = final_CD4or8 .* final_TEX/100; # the proportion of those that are CD4/CD8 * the proportion of those that are PD-1+TIM-3+

    # final amounts of each T cell subtype 
    final_num = init_T * fold_exp_comp .* final_prop; # [naive/activated/memory/effector/TPEX/TEX, 7 stimulus amounts, CD4/CD8, healthy/patient, 10 samples]




    # re-form cell count data 

    # sum over CD4 and CD8
    init_num_total = init_num[:,1,:,:] + init_num[:,2,:,:]; # [naive/activated/memory/effector/TPEX/TEX, healthy/patient, 10 samples]
    init_prop_total = init_prop[:,1,:,:] + init_prop[:,2,:,:];
    final_num_total = final_num[:,:,1,:,:] + final_num[:,:,2,:,:]; # [naive/activated/memory/effector/TPEX/TEX, 7 stimulus amounts, healthy/patient, 10 samples]
    final_prop_total = final_prop[:,:,1,:,:] + final_prop[:,:,2,:,:];

    # get means and standard deviations over all samples
    mean_init_num = zeros(6,2); mean_init_num[:,1] = mean(init_num_total[:,1,1:8], dims=2)[:,1]; mean_init_num[:,2] = mean(init_num_total[:,2,:], dims=2)[:,1]; # [naive/activated/memory/effector/TPEX/TEX, healthy/patient]
    std_init_num = zeros(6,2); std_init_num[:,1] = std(init_num_total[:,1,1:8], dims=2)[:,1]; std_init_num[:,2] = std(init_num_total[:,2,:], dims=2)[:,1]; # [naive/activated/memory/effector/TPEX/TEX, healthy/patient]
    mean_init_prop = zeros(6,2); mean_init_prop[:,1] = mean(init_prop_total[:,1,1:8], dims=2)[:,1]; mean_init_prop[:,2] = mean(init_prop_total[:,2,:], dims=2)[:,1]; # [naive/activated/memory/effector/TPEX/TEX, healthy/patient]
    std_init_prop = zeros(6,2); std_init_prop[:,1] = std(init_prop_total[:,1,1:8], dims=2)[:,1]; std_init_prop[:,2] = std(init_prop_total[:,2,:], dims=2)[:,1]; # [naive/activated/memory/effector/TPEX/TEX, healthy/patient]

    mean_final_num = zeros(6,7,2); std_final_num = zeros(6,7,2); # [naive/activated/memory/effector/TPEX/TEX, 7 stimulus amounts, healthy/patient]
    mean_final_prop = zeros(6,7,2); std_final_prop = zeros(6,7,2); 
    mean_final_totnum = zeros(7,2); std_final_totnum = zeros(7,2);  # total number of cells: [7 stimulus amounts, healthy/patient]
    for stim_amount = 1:7 
        for donor = 1:2 
            totnum = []; # initialise total number
            for cell_type = 1:6 # for each entry in the final mean/std arrays
                non_NaN_num = []; # non-NaN entries for the number of cells 
                non_NaN_prop = []; # non-NaN entries for the proportion of cells
                for sample = 1:10 # average or take standard deviation of non-NaN entries 
                    if !isnan(final_num_total[cell_type,stim_amount,donor,sample]) # if entry is not NaN, store it
                        append!(non_NaN_num, final_num_total[cell_type,stim_amount,donor,sample])
                        append!(non_NaN_prop, final_prop_total[cell_type,stim_amount,donor,sample])
                    end
                end
                mean_final_num[cell_type,stim_amount,donor] = mean(non_NaN_num); # compute mean and standard deviation of non-NaN entries
                std_final_num[cell_type,stim_amount,donor] = std(non_NaN_num);
                mean_final_prop[cell_type,stim_amount,donor] = mean(non_NaN_prop); 
                std_final_prop[cell_type,stim_amount,donor] = std(non_NaN_prop);
                if cell_type == 1
                    totnum = non_NaN_num;
                elseif cell_type == 3 || cell_type == 4 || cell_type == 6
                    totnum += non_NaN_num
                end
            end
            mean_final_totnum[stim_amount,donor] = mean(totnum); # compute mean and standard deviation of non-NaN entries
            std_final_totnum[stim_amount,donor] = std(totnum);
        end
    end



    if plot_data # plot cell count data 
        col_palette = palette(:viridis,7); # colour palette for plotting
        cell_names = ["Naive","Activated","Memory","Effector","TPEX","TEX"]; # cell names for plot titles 
        cell_cols = [col_nai,col_act,col_mem,col_eff,col_pro,col_exh]; # colours for each cell type 
        donor_names = ["healthy","patient"]; # donor type for plot titles
        for donor = 1:2 # for healthy and patient donor samples 
            # plot initial proportion of each cell type
            initprop_plt = CairoMakie.Figure(fontsize=25); ax = CairoMakie.Axis(initprop_plt[1, 1], xticks = (1:6, cell_names), xticklabelrotation=45.0, width=450, height=250, xlabel="Cell type", ylabel="Proportion of total cells", title="$(donor_names[donor])");
            CairoMakie.ylims!(0,maximum(mean_init_prop+std_init_prop)*1.05)
            if donor == 2; 
                ax.ylabel="";
                ax.yticklabelsvisible = false; 
            end
            CairoMakie.barplot!(mean_init_prop[:,donor], color=cell_cols)
            CairoMakie.errorbars!(1:6, mean_init_prop[:,donor], std_init_prop[:,donor], color=:black, whiskerwidth=10)
            display(initprop_plt)

            # plot initial numbers of each cell type
            init_plt = CairoMakie.Figure(fontsize=25); ax = CairoMakie.Axis(init_plt[1, 1], xticks = (1:6, cell_names), xticklabelrotation=45.0, width=450, height=250, xlabel="Cell type", ylabel="Number of cells", title="$(donor_names[donor])");
            CairoMakie.ylims!(0,maximum(mean_init_num+std_init_num)*1.05)
            if donor == 2; 
                ax.ylabel="";
                ax.yticklabelsvisible = false; 
            end
            CairoMakie.barplot!(mean_init_num[:,donor], color=cell_cols)
            CairoMakie.errorbars!(1:6, mean_init_num[:,donor], std_init_num[:,donor], color=:black, whiskerwidth=10)
            display(init_plt)

            # total number of cells (sum of naive, memory, effector and TEX)
            totexp_plt = CairoMakie.Figure(fontsize=25); ax = CairoMakie.Axis(totexp_plt[1, 1], width=450, height=300, xlabel="Time (days)", ylabel="Total number of cells", title="$(donor_names[donor])");
            CairoMakie.ylims!(0,maximum(mean_final_totnum+std_final_totnum)*1.05)
            if donor == 2; 
                ax.ylabel="";
                ax.yticklabelsvisible = false; 
            end
            CairoMakie.lines!([0.02,0.05,0.1,0.15,0.2,0.25,0.3], mean_final_totnum[:,donor], color=:black);
            CairoMakie.band!([0.02,0.05,0.1,0.15,0.2,0.25,0.3], mean_final_totnum[:,donor]-std_final_totnum[:,donor], mean_final_totnum[:,donor]+std_final_totnum[:,donor]; color=(:black,0.3));
            display(totexp_plt)
        end
        for cell_type = 1:6 # for each cell type
            for donor = 1:2 # for healthy/patient donor
                # plot expansion from initial to final time (stimulus value represented by colour)
                exp_plt = CairoMakie.Figure(fontsize=25); ax = CairoMakie.Axis(exp_plt[1, 1], width=450, height=300, xlabel="Time (days)", ylabel="Number of cells", title="$(cell_names[cell_type]): $(donor_names[donor])");
                for stim_val = 1:7
                    mean_init_to_final = [mean_init_num[cell_type,donor], mean_final_num[cell_type,stim_val,donor]]; # values from initial to final
                    std_init_to_final = [std_init_num[cell_type,donor], std_final_num[cell_type,stim_val,donor]]; 
                    CairoMakie.ylims!(0,maximum(mean_final_num[cell_type,:,:]+std_final_num[cell_type,:,:])*1.05)
                    if donor == 2; 
                        ax.ylabel="";
                        ax.yticklabelsvisible = false; 
                    end
                    CairoMakie.lines!([0,8], mean_init_to_final; color=col_palette[stim_val]); 
                    CairoMakie.band!([0,8], mean_init_to_final-std_init_to_final, mean_init_to_final+std_init_to_final; color=(col_palette[stim_val],0.1));
                end
                display(exp_plt)

                # plot final proportion for each stimulus value 
                perc_plt = CairoMakie.Figure(fontsize=25); ax = CairoMakie.Axis(perc_plt[1, 1], width=450, height=300, xlabel="APC-ms (mol%)", ylabel="Proportion of total cells", title="$(cell_names[cell_type]): $(donor_names[donor])");
                CairoMakie.ylims!(0,maximum(mean_final_prop[cell_type,:,:]+std_final_prop[cell_type,:,:])*1.05)
                if donor == 2; 
                    ax.ylabel="";
                    ax.yticklabelsvisible = false; 
                end
                CairoMakie.lines!([0.02,0.05,0.1,0.15,0.2,0.25,0.3], mean_final_prop[cell_type,:,donor]; color=cell_cols[cell_type]);
                CairoMakie.band!([0.02,0.05,0.1,0.15,0.2,0.25,0.3], mean_final_prop[cell_type,:,donor]-std_final_prop[cell_type,:,donor], mean_final_prop[cell_type,:,donor]+std_final_prop[cell_type,:,donor]; color=(cell_cols[cell_type],0.3));
                display(perc_plt)

                # plot final expansion for each stimulus value
                stim_plt = CairoMakie.Figure(fontsize=25); ax = CairoMakie.Axis(stim_plt[1, 1], width=450, height=300, xlabel="APC-ms (mol%)", ylabel="Number of cells", title="$(cell_names[cell_type]): $(donor_names[donor])");
                CairoMakie.ylims!(0,maximum(mean_final_num[cell_type,:,:]+std_final_num[cell_type,:,:])*1.05)
                if donor == 2; 
                    ax.ylabel="";
                    ax.yticklabelsvisible = false; 
                end
                CairoMakie.lines!([0.02,0.05,0.1,0.15,0.2,0.25,0.3], mean_final_num[cell_type,:,donor]; color=cell_cols[cell_type]);
                CairoMakie.band!([0.02,0.05,0.1,0.15,0.2,0.25,0.3], mean_final_num[cell_type,:,donor]-std_final_num[cell_type,:,donor], mean_final_num[cell_type,:,donor]+std_final_num[cell_type,:,donor]; color=(cell_cols[cell_type],0.3));
                display(stim_plt)
            end
        end
    end



    # return average and standard deviations of total cell counts at day 0 and day 8 for each stimulus value (patient-derived and healthy donors)
    return init_num_total, mean_init_num, std_init_num, final_num_total, mean_final_num, std_final_num

end

function A_and_Adash(ρ, params, funcs, derivs)
    # computes the coefficient matrix A for the ODE system, and all its partial derivatives with respect to the parameters being fitted

    A_dash = [zeros(4,4) for i in 1:length(params)]; # partial derivatives of A 
    if fit_mich_consts # if fitting Michealis constants alongside maximum rates
        if fit_prolif # also fitting proliferation rate
            r_dm_max, r_de_max, r_te_max, sρ_dm, sρ_de, sρ_te, r_p = params; # unpack current parameter values (maximum rates and Michaelis constants)
    
            A = [funcs[1](ρ,r_dm_max,r_de_max,r_te_max,sρ_dm,sρ_de,sρ_te,r_p) funcs[2](r_p) funcs[3](r_p) 0; 
                 funcs[4](ρ,r_dm_max,sρ_dm) funcs[5](ρ,r_de_max,r_te_max,sρ_de,sρ_te,r_p) funcs[6](ρ,r_dm_max,sρ_dm) 0; 
                 funcs[7](ρ,r_de_max,sρ_de) funcs[8](ρ,r_de_max,sρ_de) funcs[9](ρ,r_dm_max,r_te_max,sρ_dm,sρ_te,r_p) funcs[10](ρ,r_de_max,sρ_de); 
                 funcs[11](ρ,r_te_max,sρ_te) funcs[12](ρ,r_te_max,sρ_te) funcs[13](ρ,r_te_max,sρ_te) funcs[14](ρ,r_de_max,sρ_de)]; # coefficients in ODE system
            
            A_dash[1] = [derivs[1,1](ρ,sρ_dm) derivs[1,2] derivs[1,3] 0;
                         derivs[1,4](ρ,sρ_dm) derivs[1,5] derivs[1,6](ρ,sρ_dm) 0;
                         derivs[1,7] derivs[1,8] derivs[1,9](ρ,sρ_dm) derivs[1,10];
                         derivs[1,11] derivs[1,12] derivs[1,13] derivs[1,14]]; # partial deriv w.r.t. r_dm_max 
            A_dash[2] = [derivs[2,1](ρ,sρ_de) derivs[2,2] derivs[2,3] 0;
                         derivs[2,4] derivs[2,5](ρ,sρ_de) derivs[2,6] 0;
                         derivs[2,7](ρ,sρ_de) derivs[2,8](ρ,sρ_de) derivs[2,9] derivs[2,10](ρ,sρ_de);
                         derivs[2,11] derivs[2,12] derivs[2,13] derivs[2,14](ρ,sρ_de)]; # partial deriv w.r.t. r_de_max 
            A_dash[3] = [derivs[3,1](ρ,sρ_te) derivs[3,2] derivs[3,3] 0;
                         derivs[3,4] derivs[3,5](ρ,sρ_te) derivs[3,6] 0;
                         derivs[3,7] derivs[3,8] derivs[3,9](ρ,sρ_te) derivs[3,10];
                         derivs[3,11](ρ,sρ_te) derivs[3,12](ρ,sρ_te) derivs[3,13](ρ,sρ_te) derivs[3,14]]; # partial deriv w.r.t. r_te_max 
            A_dash[4] = [derivs[4,1](ρ,r_dm_max,sρ_dm) derivs[4,2] derivs[4,3] 0;
                         derivs[4,4](ρ,r_dm_max,sρ_dm) derivs[4,5] derivs[4,6](ρ,r_dm_max,sρ_dm) 0;
                         derivs[4,7] derivs[4,8] derivs[4,9](ρ,r_dm_max,sρ_dm) derivs[4,10];
                         derivs[4,11] derivs[4,12] derivs[4,13] derivs[4,14]]; # partial deriv w.r.t. sρ_dm
            A_dash[5] = [derivs[5,1](ρ,r_de_max,sρ_de) derivs[5,2] derivs[5,3] 0;
                         derivs[5,4] derivs[5,5](ρ,r_de_max,sρ_de) derivs[5,6] 0;
                         derivs[5,7](ρ,r_de_max,sρ_de) derivs[5,8](ρ,r_de_max,sρ_de) derivs[5,9] derivs[5,10](ρ,r_de_max,sρ_de);
                         derivs[5,11] derivs[5,12] derivs[5,13] derivs[5,14](ρ,r_de_max,sρ_de)]; # partial deriv w.r.t. sρ_de
            A_dash[6] = [derivs[6,1](ρ,r_te_max,sρ_te) derivs[6,2] derivs[6,3] 0;
                         derivs[6,4] derivs[6,5](ρ,r_te_max,sρ_te) derivs[6,6] 0;
                         derivs[6,7] derivs[6,8] derivs[6,9](ρ,r_te_max,sρ_te) derivs[6,10];
                         derivs[6,11](ρ,r_te_max,sρ_te) derivs[6,12](ρ,r_te_max,sρ_te) derivs[6,13](ρ,r_te_max,sρ_te) derivs[6,14]]; # partial deriv w.r.t. sρ_te
            A_dash[7] = [derivs[7,1] derivs[7,2] derivs[7,3] 0;
                         derivs[7,4] derivs[7,5] derivs[7,6] 0;
                         derivs[7,7] derivs[7,8] derivs[7,9] derivs[7,10];
                         derivs[7,11] derivs[7,12] derivs[7,13] derivs[7,14]]; # partial deriv w.r.t. r_p

        else # fitting all parameters for diff and exh rate functions
            r_dm_max, r_de_max, r_te_max, sρ_dm, sρ_de, sρ_te = params; # unpack current parameter values (maximum rates and Michaelis constants)
    
            A = [funcs[1](ρ,r_dm_max,r_de_max,r_te_max,sρ_dm,sρ_de,sρ_te) funcs[2] funcs[3] 0; 
                 funcs[4](ρ,r_dm_max,sρ_dm) funcs[5](ρ,r_de_max,r_te_max,sρ_de,sρ_te) funcs[6](ρ,r_dm_max,sρ_dm) 0; 
                 funcs[7](ρ,r_de_max,sρ_de) funcs[8](ρ,r_de_max,sρ_de) funcs[9](ρ,r_dm_max,r_te_max,sρ_dm,sρ_te) funcs[10](ρ,r_de_max,sρ_de); 
                 funcs[11](ρ,r_te_max,sρ_te) funcs[12](ρ,r_te_max,sρ_te) funcs[13](ρ,r_te_max,sρ_te) funcs[14](ρ,r_de_max,sρ_de)]; # coefficients in ODE system
            
            A_dash[1] = [derivs[1,1](ρ,sρ_dm) derivs[1,2] derivs[1,3] 0;
                         derivs[1,4](ρ,sρ_dm) derivs[1,5] derivs[1,6](ρ,sρ_dm) 0;
                         derivs[1,7] derivs[1,8] derivs[1,9](ρ,sρ_dm) derivs[1,10];
                         derivs[1,11] derivs[1,12] derivs[1,13] derivs[1,14]]; # partial deriv w.r.t. r_dm_max 
            A_dash[2] = [derivs[2,1](ρ,sρ_de) derivs[2,2] derivs[2,3] 0;
                         derivs[2,4] derivs[2,5](ρ,sρ_de) derivs[2,6] 0;
                         derivs[2,7](ρ,sρ_de) derivs[2,8](ρ,sρ_de) derivs[2,9] derivs[2,10](ρ,sρ_de);
                         derivs[2,11] derivs[2,12] derivs[2,13] derivs[2,14](ρ,sρ_de)]; # partial deriv w.r.t. r_de_max 
            A_dash[3] = [derivs[3,1](ρ,sρ_te) derivs[3,2] derivs[3,3] 0;
                         derivs[3,4] derivs[3,5](ρ,sρ_te) derivs[3,6] 0;
                         derivs[3,7] derivs[3,8] derivs[3,9](ρ,sρ_te) derivs[3,10];
                         derivs[3,11](ρ,sρ_te) derivs[3,12](ρ,sρ_te) derivs[3,13](ρ,sρ_te) derivs[3,14]]; # partial deriv w.r.t. r_te_max 
            A_dash[4] = [derivs[4,1](ρ,r_dm_max,sρ_dm) derivs[4,2] derivs[4,3] 0;
                         derivs[4,4](ρ,r_dm_max,sρ_dm) derivs[4,5] derivs[4,6](ρ,r_dm_max,sρ_dm) 0;
                         derivs[4,7] derivs[4,8] derivs[4,9](ρ,r_dm_max,sρ_dm) derivs[4,10];
                         derivs[4,11] derivs[4,12] derivs[4,13] derivs[4,14]]; # partial deriv w.r.t. sρ_dm
            A_dash[5] = [derivs[5,1](ρ,r_de_max,sρ_de) derivs[5,2] derivs[5,3] 0;
                         derivs[5,4] derivs[5,5](ρ,r_de_max,sρ_de) derivs[5,6] 0;
                         derivs[5,7](ρ,r_de_max,sρ_de) derivs[5,8](ρ,r_de_max,sρ_de) derivs[5,9] derivs[5,10](ρ,r_de_max,sρ_de);
                         derivs[5,11] derivs[5,12] derivs[5,13] derivs[5,14](ρ,r_de_max,sρ_de)]; # partial deriv w.r.t. sρ_de
            A_dash[6] = [derivs[6,1](ρ,r_te_max,sρ_te) derivs[6,2] derivs[6,3] 0;
                         derivs[6,4] derivs[6,5](ρ,r_te_max,sρ_te) derivs[6,6] 0;
                         derivs[6,7] derivs[6,8] derivs[6,9](ρ,r_te_max,sρ_te) derivs[6,10];
                         derivs[6,11](ρ,r_te_max,sρ_te) derivs[6,12](ρ,r_te_max,sρ_te) derivs[6,13](ρ,r_te_max,sρ_te) derivs[6,14]]; # partial deriv w.r.t. sρ_te
        end

    else # not fitting Michaelis constants
        if fit_prolif # fitting maximum rates for diff and exh, and proliferation rate
            r_dm_max, r_de_max, r_te_max, r_p = params; # only unpack maximum rates
    
            A = [funcs[1](ρ,r_dm_max,r_de_max,r_te_max,r_p) funcs[2](r_p) funcs[3](r_p) 0; 
                 funcs[4](ρ,r_dm_max) funcs[5](ρ,r_de_max,r_te_max,r_p) funcs[6](ρ,r_dm_max) 0; 
                 funcs[7](ρ,r_de_max) funcs[8](ρ,r_de_max) funcs[9](ρ,r_dm_max,r_te_max,r_p) funcs[10](ρ,r_de_max); 
                 funcs[11](ρ,r_te_max) funcs[12](ρ,r_te_max) funcs[13](ρ,r_te_max) funcs[14](ρ,r_de_max)]; # coefficients in ODE system
            
            A_dash[1] = [derivs[1,1](ρ) derivs[1,2] derivs[1,3] 0;
                         derivs[1,4](ρ) derivs[1,5] derivs[1,6](ρ) 0;
                         derivs[1,7] derivs[1,8] derivs[1,9](ρ) derivs[1,10];
                         derivs[1,11] derivs[1,12] derivs[1,13] derivs[1,14]]; # partial deriv w.r.t. r_dm_max 
            A_dash[2] = [derivs[2,1](ρ) derivs[2,2] derivs[2,3] 0;
                         derivs[2,4] derivs[2,5](ρ) derivs[2,6] 0;
                         derivs[2,7](ρ) derivs[2,8](ρ) derivs[2,9] derivs[2,10](ρ);
                         derivs[2,11] derivs[2,12] derivs[2,13] derivs[2,14](ρ)]; # partial deriv w.r.t. r_de_max 
            A_dash[3] = [derivs[3,1](ρ) derivs[3,2] derivs[3,3] 0;
                         derivs[3,4] derivs[3,5](ρ) derivs[3,6] 0;
                         derivs[3,7] derivs[3,8] derivs[3,9](ρ) derivs[3,10];
                         derivs[3,11](ρ) derivs[3,12](ρ) derivs[3,13](ρ) derivs[3,14]]; # partial deriv w.r.t. r_te_max 
            A_dash[4] = [derivs[4,1] derivs[4,2] derivs[4,3] 0;
                         derivs[4,4] derivs[4,5] derivs[4,6] 0;
                         derivs[4,7] derivs[4,8] derivs[4,9] derivs[4,10];
                         derivs[4,11] derivs[4,12] derivs[4,13] derivs[4,14]]; # partial deriv w.r.t. r_p 

        else # fitting only maximum rates for diff and exh
            r_dm_max, r_de_max, r_te_max = params; # only unpack maximum rates
    
            A = [funcs[1](ρ,r_dm_max,r_de_max,r_te_max) funcs[2] funcs[3] 0; 
                 funcs[4](ρ,r_dm_max) funcs[5](ρ,r_de_max,r_te_max) funcs[6](ρ,r_dm_max) 0; 
                 funcs[7](ρ,r_de_max) funcs[8](ρ,r_de_max) funcs[9](ρ,r_dm_max,r_te_max) funcs[10](ρ,r_de_max); 
                 funcs[11](ρ,r_te_max) funcs[12](ρ,r_te_max) funcs[13](ρ,r_te_max) funcs[14](ρ,r_de_max)]; # coefficients in ODE system
            
            A_dash[1] = [derivs[1,1](ρ) derivs[1,2] derivs[1,3] 0;
                         derivs[1,4](ρ) derivs[1,5] derivs[1,6](ρ) 0;
                         derivs[1,7] derivs[1,8] derivs[1,9](ρ) derivs[1,10];
                         derivs[1,11] derivs[1,12] derivs[1,13] derivs[1,14]]; # partial deriv w.r.t. r_dm_max 
            A_dash[2] = [derivs[2,1](ρ) derivs[2,2] derivs[2,3] 0;
                         derivs[2,4] derivs[2,5](ρ) derivs[2,6] 0;
                         derivs[2,7](ρ) derivs[2,8](ρ) derivs[2,9] derivs[2,10](ρ);
                         derivs[2,11] derivs[2,12] derivs[2,13] derivs[2,14](ρ)]; # partial deriv w.r.t. r_de_max 
            A_dash[3] = [derivs[3,1](ρ) derivs[3,2] derivs[3,3] 0;
                         derivs[3,4] derivs[3,5](ρ) derivs[3,6] 0;
                         derivs[3,7] derivs[3,8] derivs[3,9](ρ) derivs[3,10];
                         derivs[3,11](ρ) derivs[3,12](ρ) derivs[3,13](ρ) derivs[3,14]]; # partial deriv w.r.t. r_te_max 
        end
    end

    return A, A_dash;
end

function sol_and_grad(s, u0, T, ρ, params, funcs, derivs, scaler)
    # computes the solution to the ODE system at time T, along with the derivatives of the solution with respect to each parameter being fitted

    A, A_dash = A_and_Adash(s*ρ, params, funcs, derivs); # get the matrix A and its derivatives w.r.t. parameters in params 

    uT = exp(A*T)*u0/scaler; # scaled down solution at time T:  [naive, memory, effector, exhausted] 
    
    ∂u = [zeros(length(u0)) for i = eachindex(params)]; # gradient of solution for each parameter in params 
    
    n = length(u0); # number of cell types
    M = zeros(2*n, 2*n); # block matrix for exponentiation
    M[1:n, 1:n] = A*T; # fill top left and bottom right corners of block matrix with A*T
    M[n+1:2*n, n+1:2*n] = A*T;
    for k = eachindex(params) # for each parameter
        M[1:n, n+1:2*n] = A_dash[k]*T; # fill top right of block matrix with the partial derivative of A w.r.t. the kth parameter (multiplied by T)
        ∂u[k] = exp(M)[1:n, n+1:2*n]*u0/scaler; # the kth partial derivative of the solution is the upper right corner of e^M multiplied by u0
    end

    return uT, ∂u
end 

function loss_and_grad!(grad, params, data, sol, sol_grad, stim_vals, donors)
    # computes the loss function used for optimisation, along with the gradient of the loss function with respect to the parameters being fitted

    L = 0; # initialise loss value for cell count
    curr_grad = zeros(size(grad)); # initialise current gradient of the loss function
    for s = stim_vals # for each stimulus value
        stim_ind = findfirst(x->x==s,[0.8,2,4,6,8,10,12]); # index for current stimulus ratio 
        for donor_n = donors # for each donor ID 
            u0 = ICs[donor_n]; # initial condition [naive, memory, effector, exhausted]

            u_i = sol(s, u0, params); # solution to ODE system for current set of parameters/IC 
            data_i = data[:,stim_ind,donor_n]; # current data

            err_i = u_i - data_i; # error to data for current set of parameters

            L += dot(err_i, err_i); # add the sum of squared error to the current loss

            ∂u_i = sol_grad(s, u0, params); # gradient of the solution (∂u/∂pₖ) and gradients for current set of parameters
            for k = eachindex(grad) # for each parameter
                curr_grad[k] += dot(∂u_i[k], err_i); # add to the sum in the equation for gradient of the loss
            end
        end
    end
    grad .= 2*curr_grad/(length(donors)*length(stim_vals)); # update gradient of the loss function
    
    return L/(length(donors)*length(stim_vals)); # return the average sum of squared errors as the loss
end

function plot_pathway(model_n, c1a, c1b, c1c, c1d, c1e, c1f, c1g, c2a, c2b, c2c, c3a, c3b, c3c, c3d, c3e)
    # plots the differentiation and exhaustion pathway for the current set of model choices (models with naive, memory, effector, exhausted cells - excluding choices 1f and 1i)
    
    choice_names = ["1a","1b","1c","1d","1e","1f","1g","2a","2b","2c","3a","3b","3c","3d","3e"]; # names for each model choice

    col_diff = RGB(0/255,0/255,0/255); # colour for state transitions/differentiation and exhaustion arrows
    col_prolif = RGB(33/255,69/255,166/255); # colour for proliferation event arrows
    col_death = RGB(145/255,33/255,33/255); # colour for death event arrows
    arrow_width = 5; # linewidth for arrows

    cell_pos = [-1 0; 0 0; 1 0; 0.5 -0.5]; # positions of each cell in the graph [cell type, x/y]  (origin on memory cell)

    state_graph = CairoMakie.Figure(fontsize=25, size=(600,350)); # initialise graph
    ax = CairoMakie.Axis(state_graph[1, 1], limits=(minimum(cell_pos[:,1])-0.25,maximum(cell_pos[:,1])+0.26,minimum(cell_pos[:,2])-0.26,maximum(cell_pos[:,2])+0.4), aspect=DataAspect(), title="Model $(model_n):\n$(join(["$(choice_names[i])" for i in findall(x->x==true,[c1a,c1b,c1c,c1d,c1e,c1f,c1g,c2a,c2b,c2c,c3a,c3b,c3c,c3d,c3e])],", ")) on, the rest off"); # define axis with title
    hidespines!(ax); # remove axis grid and other background elements
    hidedecorations!(ax);
    set_theme!(fonts = (; regular = "Times New Roman", bold = "Times New Roman Bold")); # set fonts

    CairoMakie.scatter!(cell_pos[:,1], cell_pos[:,2]; markersize=100, color=[col_nai,col_mem,col_eff,col_exh], strokecolor=:black, strokewidth=3)
    CairoMakie.scatter!(cell_pos[:,1]+(rand(4).-0.5)/20, cell_pos[:,2]+(rand(4).-0.5)/20; markersize=60, color=[col_nai,col_mem,col_eff,col_exh]/1.5)

    if use_arrows2d # use arrows2d function
        # all state transitions: differentiation and exhaustion 
        if (!c1a && !c1c) || c1b; CairoMakie.lines!([cell_pos[1,1]+0.2,cell_pos[2,1]-0.22], [cell_pos[1,2],cell_pos[2,2]], color=col_diff, linewidth=arrow_width); CairoMakie.arrows2d!([cell_pos[2,1]-0.26], [cell_pos[2,2]], [0.1], [0]; tipwidth=arrow_width*4, tiplength=arrow_width*3, shaftwidth=0); end # line and arrowhead for naive -> memory
        if c1a || c1b || c1c; CairoMakie.lines!(cell_pos[1,1]+0.14:0.01:cell_pos[3,1]-0.15, -0.11*(cell_pos[1,1]+0.14:0.01:cell_pos[3,1]-0.15).^2 .+ 0.22, color=col_diff, linewidth=arrow_width); CairoMakie.arrows2d!([cell_pos[3,1]-0.19], [-0.11*(cell_pos[3,1]-0.19)^2+0.225], [0.1], [-0.04]; tipwidth=arrow_width*4, tiplength=arrow_width*3, shaftwidth=0); end # naive -> effector
        if c2a; CairoMakie.lines!(cell_pos[1,1]+0.16:0.01:cell_pos[4,1]-0.22, 0.29*((cell_pos[1,1]+0.16:0.01:cell_pos[4,1]-0.22).-0.3).^2 .- 0.5, color=col_diff, linewidth=arrow_width); CairoMakie.arrows2d!([cell_pos[4,1]-0.26], [0.29*(cell_pos[4,1]-0.54)^2-0.5], [0.1], [0]; tipwidth=arrow_width*4, tiplength=arrow_width*3, shaftwidth=0); end # naive -> exhausted
        if (!c1a || c1c) && !c1d; CairoMakie.lines!([cell_pos[2,1]+0.2,cell_pos[3,1]-0.22], [cell_pos[2,2]+0.03,cell_pos[3,2]+0.03], color=col_diff, linewidth=arrow_width); CairoMakie.arrows2d!([cell_pos[3,1]-0.26], [cell_pos[3,2]+0.03], [0.1], [0]; tipwidth=arrow_width*4, tiplength=arrow_width*3, shaftwidth=0); end # memory -> effector
        if c2b; CairoMakie.lines!([cell_pos[2,1]+0.14,cell_pos[4,1]-0.15], [cell_pos[2,2]-0.14,cell_pos[4,2]+0.15], color=col_diff, linewidth=arrow_width); CairoMakie.arrows2d!([cell_pos[4,1]-0.18], [cell_pos[4,2]+0.18], [0.1], [-0.1]; tipwidth=arrow_width*4, tiplength=arrow_width*3, shaftwidth=0); end # memory -> exhausted
        if (c1a || c1c) && !c1d; CairoMakie.lines!([cell_pos[2,1]+0.22,cell_pos[3,1]-0.2], [cell_pos[2,2]-0.03,cell_pos[3,2]-0.03], color=col_diff, linewidth=arrow_width); CairoMakie.arrows2d!([cell_pos[2,1]+0.26], [cell_pos[2,2]-0.03], [-0.1], [0]; tipwidth=arrow_width*4, tiplength=arrow_width*3, shaftwidth=0); end # effector -> memory
        if c2c; CairoMakie.lines!([cell_pos[3,1]-0.125,cell_pos[4,1]+0.165], [cell_pos[3,2]-0.165,cell_pos[4,2]+0.125], color=col_diff, linewidth=arrow_width); CairoMakie.arrows2d!([cell_pos[4,1]+0.2], [cell_pos[4,2]+0.16], [-0.1], [-0.1]; tipwidth=arrow_width*4, tiplength=arrow_width*3, shaftwidth=0); end # effector -> exhausted
        if c1g; CairoMakie.lines!([cell_pos[3,1]-0.165,cell_pos[4,1]+0.125], [cell_pos[3,2]-0.125,cell_pos[4,2]+0.165], color=col_diff, linewidth=arrow_width); CairoMakie.arrows2d!([cell_pos[3,1]-0.2], [cell_pos[3,2]-0.16], [0.1], [0.1]; tipwidth=arrow_width*4, tiplength=arrow_width*3, shaftwidth=0); end # exhausted -> effector

        # all proliferation events       used to have  linestyle=(:dash,1.5)
        if c3c; CairoMakie.arc!(cell_pos[1,:]+[-0.1,0.2], 0.12, 1.3*pi, -0.06*pi, color=col_prolif, linewidth=arrow_width); CairoMakie.arrows2d!([cell_pos[1,1]-0.2], [cell_pos[1,2]+0.13], [0.1], [-0.08]; tipwidth=arrow_width*4, tiplength=arrow_width*3, shaftwidth=0, color=col_prolif); end # symmetric naive proliferation
        if !c3a; CairoMakie.arc!(cell_pos[2,:]+[-0.05,-0.22], 0.12, -1.25*pi, 0.05*pi, color=col_prolif, linewidth=arrow_width); CairoMakie.arrows2d!([cell_pos[2,1]+0.065], [cell_pos[2,2]-0.23], [-0.02], [0.1]; tipwidth=arrow_width*4, tiplength=arrow_width*3, shaftwidth=0, color=col_prolif); end # symmetric memory proliferation
        if c3a; CairoMakie.lines!(cell_pos[1,1]+0.2:0.01:cell_pos[2,1]-0.18, 0.3*((cell_pos[1,1]+0.2:0.01:cell_pos[2,1]-0.18).+0.5).^2 .- 0.09, color=col_prolif, linewidth=arrow_width); CairoMakie.arrows2d!([cell_pos[1,1]+0.24], [0.3*(cell_pos[1,1]+0.74)^2-0.09], [-0.1], [0.02]; tipwidth=arrow_width*4, tiplength=arrow_width*3, shaftwidth=0, color=col_prolif); end # asymmetric memory proliferation
        if !c3a && c3b; CairoMakie.arc!(cell_pos[3,:]+[0.13,0.2], 0.12, 1.0*pi, -0.4*pi, color=col_prolif, linewidth=arrow_width); CairoMakie.arrows2d!([cell_pos[3,1]+0.01], [cell_pos[3,2]+0.24], [0], [-0.1]; tipwidth=arrow_width*4, tiplength=arrow_width*3, shaftwidth=0, color=col_prolif); end # symmetric effector proliferation
        if c3a && c3b; CairoMakie.lines!(cell_pos[1,1]+0.1:0.01:cell_pos[3,1]-0.08, -0.17*(cell_pos[1,1]+0.1:0.01:cell_pos[3,1]-0.08).^2 .+ 0.33, color=col_prolif, linewidth=arrow_width); CairoMakie.arrows2d!([cell_pos[1,1]+0.13], [-0.17*(cell_pos[1,1]+0.13)^2+0.33], [-0.1], [-0.07]; tipwidth=arrow_width*4, tiplength=arrow_width*3, shaftwidth=0, color=col_prolif); end # asymmetric effector proliferation

        # all death events       used to have  linestyle=(:dash,1.5)
        if c3d; CairoMakie.lines!([cell_pos[3,1]+0.13,cell_pos[3,1]+0.25], [cell_pos[3,2]-0.13,cell_pos[3,2]-0.25], color=col_death, linewidth=arrow_width); CairoMakie.arrows2d!([cell_pos[3,1]+0.215], [cell_pos[3,2]-0.215], [0.1], [-0.1]; tipwidth=arrow_width*4, tiplength=arrow_width*3, shaftwidth=0, color=col_death); end # effector death
        if c3e; CairoMakie.lines!([cell_pos[4,1]+0.13,cell_pos[4,1]+0.25], [cell_pos[4,2]-0.13,cell_pos[4,2]-0.25], color=col_death, linewidth=arrow_width); CairoMakie.arrows2d!([cell_pos[4,1]+0.215], [cell_pos[4,2]-0.215], [0.1], [-0.1]; tipwidth=arrow_width*4, tiplength=arrow_width*3, shaftwidth=0, color=col_death); end # exhausted death
    else # use arrows function
        # all state transitions: differentiation and exhaustion 
        if (!c1a && !c1c) || c1b; CairoMakie.lines!([cell_pos[1,1]+0.2,cell_pos[2,1]-0.22], [cell_pos[1,2],cell_pos[2,2]], color=col_diff, linewidth=arrow_width); CairoMakie.arrows!([cell_pos[2,1]-0.24], [cell_pos[2,2]], [0.01], [0]; arrowsize = arrow_width*4); end # line and arrowhead for naive -> memory
        if c1a || c1b || c1c; CairoMakie.lines!(cell_pos[1,1]+0.14:0.01:cell_pos[3,1]-0.15, -0.11*(cell_pos[1,1]+0.14:0.01:cell_pos[3,1]-0.15).^2 .+ 0.22, color=col_diff, linewidth=arrow_width); CairoMakie.arrows!([cell_pos[3,1]-0.175], [-0.11*(cell_pos[3,1]-0.175)^2+0.22], [0.01], [-0.004]; arrowsize = arrow_width*4); end # naive -> effector
        if c2a; CairoMakie.lines!(cell_pos[1,1]+0.16:0.01:cell_pos[4,1]-0.22, 0.29*((cell_pos[1,1]+0.16:0.01:cell_pos[4,1]-0.22).-0.3).^2 .- 0.5, color=col_diff, linewidth=arrow_width); CairoMakie.arrows!([cell_pos[4,1]-0.24], [0.29*(cell_pos[4,1]-0.52)^2-0.5], [0.01], [0]; arrowsize = arrow_width*4); end # naive -> exhausted
        if (!c1a || c1c) && !c1d; CairoMakie.lines!([cell_pos[2,1]+0.2,cell_pos[3,1]-0.22], [cell_pos[2,2]+0.03,cell_pos[3,2]+0.03], color=col_diff, linewidth=arrow_width); CairoMakie.arrows!([cell_pos[3,1]-0.24], [cell_pos[3,2]+0.03], [0.01], [0]; arrowsize = arrow_width*4); end # memory -> effector
        if c2b; CairoMakie.lines!([cell_pos[2,1]+0.14,cell_pos[4,1]-0.15], [cell_pos[2,2]-0.14,cell_pos[4,2]+0.15], color=col_diff, linewidth=arrow_width); CairoMakie.arrows!([cell_pos[4,1]-0.17], [cell_pos[4,2]+0.17], [0.01], [-0.01]; arrowsize = arrow_width*4); end # memory -> exhausted
        if (c1a || c1c) && !c1d; CairoMakie.lines!([cell_pos[2,1]+0.22,cell_pos[3,1]-0.2], [cell_pos[2,2]-0.03,cell_pos[3,2]-0.03], color=col_diff, linewidth=arrow_width); CairoMakie.arrows!([cell_pos[2,1]+0.24], [cell_pos[2,2]-0.03], [-0.01], [0]; arrowsize = arrow_width*4); end # effector -> memory
        if c2c; CairoMakie.lines!([cell_pos[3,1]-0.125,cell_pos[4,1]+0.165], [cell_pos[3,2]-0.165,cell_pos[4,2]+0.125], color=col_diff, linewidth=arrow_width); CairoMakie.arrows!([cell_pos[4,1]+0.19], [cell_pos[4,2]+0.15], [-0.01], [-0.01]; arrowsize = arrow_width*4); end # effector -> exhausted
        if c1g; CairoMakie.lines!([cell_pos[3,1]-0.165,cell_pos[4,1]+0.125], [cell_pos[3,2]-0.125,cell_pos[4,2]+0.165], color=col_diff, linewidth=arrow_width); CairoMakie.arrows!([cell_pos[3,1]-0.19], [cell_pos[3,2]-0.15], [0.01], [0.01]; arrowsize = arrow_width*4); end # exhausted -> effector

        # all proliferation events       used to have  linestyle=(:dash,1.5)
        if c3c; CairoMakie.arc!(cell_pos[1,:]+[-0.1,0.2], 0.12, 1.3*pi, -0.06*pi, color=col_prolif, linewidth=arrow_width); CairoMakie.arrows!([cell_pos[1,1]-0.19], [cell_pos[1,2]+0.12], [0.01], [-0.008]; arrowsize = arrow_width*4, color=col_prolif); end # symmetric naive proliferation
        if !c3a; CairoMakie.arc!(cell_pos[2,:]+[-0.05,-0.22], 0.12, -1.25*pi, 0.05*pi, color=col_prolif, linewidth=arrow_width); CairoMakie.arrows!([cell_pos[2,1]+0.065], [cell_pos[2,2]-0.21], [-0.0002], [0.001]; arrowsize = arrow_width*4, color=col_prolif); end # symmetric memory proliferation
        if c3a; CairoMakie.lines!(cell_pos[1,1]+0.2:0.01:cell_pos[2,1]-0.18, 0.3*((cell_pos[1,1]+0.2:0.01:cell_pos[2,1]-0.18).+0.5).^2 .- 0.09, color=col_prolif, linewidth=arrow_width); CairoMakie.arrows!([cell_pos[1,1]+0.23], [0.3*(cell_pos[1,1]+0.73)^2-0.09], [-0.01], [0.002]; arrowsize = arrow_width*4, color=col_prolif); end # asymmetric memory proliferation
        if !c3a && c3b; CairoMakie.arc!(cell_pos[3,:]+[0.13,0.2], 0.12, 1.0*pi, -0.4*pi, color=col_prolif, linewidth=arrow_width); CairoMakie.arrows!([cell_pos[3,1]+0.01], [cell_pos[3,2]+0.22], [0.0], [-0.001]; arrowsize = arrow_width*4, color=col_prolif); end # symmetric effector proliferation
        if c3a && c3b; CairoMakie.lines!(cell_pos[1,1]+0.1:0.01:cell_pos[3,1]-0.08, -0.17*(cell_pos[1,1]+0.1:0.01:cell_pos[3,1]-0.08).^2 .+ 0.33, color=col_prolif, linewidth=arrow_width); CairoMakie.arrows!([cell_pos[1,1]+0.12], [-0.17*(cell_pos[1,1]+0.12)^2+0.33], [-0.01], [-0.007]; arrowsize = arrow_width*4, color=col_prolif); end # asymmetric effector proliferation

        # all death events       used to have  linestyle=(:dash,1.5)
        if c3d; CairoMakie.lines!([cell_pos[3,1]+0.13,cell_pos[3,1]+0.25], [cell_pos[3,2]-0.13,cell_pos[3,2]-0.25], color=col_death, linewidth=arrow_width); CairoMakie.arrows!([cell_pos[3,1]+0.22], [cell_pos[3,2]-0.22], [0.01], [-0.01]; arrowsize = arrow_width*4, color=col_death); end # effector death
        if c3e; CairoMakie.lines!([cell_pos[4,1]+0.13,cell_pos[4,1]+0.25], [cell_pos[4,2]-0.13,cell_pos[4,2]-0.25], color=col_death, linewidth=arrow_width); CairoMakie.arrows!([cell_pos[4,1]+0.22], [cell_pos[4,2]-0.22], [0.01], [-0.01]; arrowsize = arrow_width*4, color=col_death); end # exhausted death
    end

    # legend 
    CairoMakie.scatter!([-2], [-2]; markersize=25, color=col_nai, strokecolor=:black, strokewidth=2, label="Naive")
    CairoMakie.scatter!([-2], [-2]; markersize=25, color=col_mem, strokecolor=:black, strokewidth=2, label="Memory")
    CairoMakie.scatter!([-2], [-2]; markersize=25, color=col_eff, strokecolor=:black, strokewidth=2, label="Effector")
    CairoMakie.scatter!([-2], [-2]; markersize=25, color=col_exh, strokecolor=:black, strokewidth=2, label="Exhausted")
    CairoMakie.lines!([-2], [-2], color=col_diff, linewidth=arrow_width*0.75, label="Transition");
    CairoMakie.lines!([-2], [-2], color=col_prolif, linewidth=arrow_width*0.75, label="Proliferation"); # used to have  linestyle=(:dot,1.2)
    CairoMakie.lines!([-2], [-2], color=col_death, linewidth=arrow_width*0.75, label="Death"); # used to have  linestyle=(:dot,1.2)
    CairoMakie.axislegend(position=(0,-0.175), framevisible=false, nbanks=2, labelsize=18);

    display(state_graph); # view graph
    return state_graph
end