#here is the OMI analysis file

using Distributions, FFTW, StatsBase, Plots, ImageView, Images, OffsetArrays, CoordinateTransformations, Dierckx, QuartzImageIO, Statistics

##### load image #####

"""
    ft_filter(filename:string)

Load the image file and convert it to fourier domain to filter
"""

function ft_filter(filename::String)
    #read in image
    img::Array{Float64} = load(filename) .|> Float64
    img = img./maximum(img)

    #apply mean filter
    kernel = ones(5, 5)./9 # mean filter
    img = imfilter(img, centered(kernel))

    #take fourier transform
    ft = fft(img) |> fftshift

    #remove background noise
    for i in eachindex(ft)
        if abs.(ft[i]) < 300.0
            ft[i] = 0
        end
    end

    #take inverse fourier transform
    img_processed = abs.(ifft(ft))
end


function center_img(img)
    kernel = ones(5, 5)./9
    dimg = imfilter(img,centered(kernel))
    gy, gx = imgradients(dimg,KernelFactors.ando5)
    angles = rad2deg.(atan.(gx ./ gy))
    d = fit_mle(Laplace, angles[:])
    theta = d.μ
    scale_factor = d.θ

    tfm = recenter(RotMatrix(deg2rad(theta)),[size(img)[1]/2,size(img)[1]/2])
    rot_img = warp(img, tfm)
    return rot_img, theta

end

function recenter_img(img,θ)
    b = axes(img)[1]
    c = []
    for i in b
       push!(c,i)
     end
     b1 = axes(img)[2]
     c1 = []
     for i in b1
        push!(c1,i)
      end
      in1 = c[Int(floor(length(c)/2))]
      in2 = c[Int(floor(length(c1)/2))]

    tfm = recenter(RotMatrix(deg2rad(-θ)),[in1,in2])
    rot_img = warp(img, tfm)
end


#get sarcomere length estimates from local maximums of columns
function get_sarclength(rot_img)

    sl_data = similar(rot_img)
    for in_col in axes(rot_img)[2]

        darray = rot_img[:,in_col]

        locmax = findlocalmaxima(darray)
        idx = []

        for i in eachindex(locmax)
            push!(idx,getindex(locmax,i)[1])
        end

        n = length(idx)

        if n > 4

            sarclength = zeros(n)
            sarclength[1] = idx[2] - idx[1]
            sarclength[n] = idx[n] - idx[n-1]
            sarclength[2:n-1] = (idx[3:n] - idx[2:n-1] +idx[2:n-1] - idx[1:n-2])/2

            spl = Spline1D(idx,sarclength)
            sl = round.(spl(Vector(idx[1]:1:idx[end])),digits=0)

            sl1 = darray
            sl1[locmax[1][1]:locmax[end][1]] = sl
        else
            sl1 = darray
        end


        for i in eachindex(sl1)
            if sl1[i] < 1.5
                sl1[i] = 0
            end
        end

    sl_data[:,in_col] = sl1
    end

return sl_data
end


function damage_kernel(test)

    freq_sl = mode(test)
    std_sl = std(test)
    weight = 1.1

    sl_datarow = size(test)[1]
    sl_datacol = size(test)[2]


    kernel = 9;
    kernel_mtrx = ones(kernel,kernel);
    kernel_mtrx2 = -1 .* kernel_mtrx;


    ctr_pxl = median(1:kernel) .|> Int64
    nbr_pxl = (ctr_pxl - 1) .|> Int64


    row_idx = range(ctr_pxl, step=kernel, stop=(sl_datarow-nbr_pxl))
    column_idx = range(ctr_pxl, step=kernel, stop=(sl_datacol-nbr_pxl))


    for col in values(column_idx)
        for row in values(row_idx)

            discrete_kernel = count(i->(i>freq_sl*weight), test[row-nbr_pxl:row+nbr_pxl, col-nbr_pxl:col+nbr_pxl]);

            if discrete_kernel > floor(kernel/2)
                test[row-nbr_pxl:row+nbr_pxl, col-nbr_pxl:col+nbr_pxl] = kernel_mtrx2;

            end
        end
    end
    return test
end


function count_dam(test, count_dmg, count_edge, count_norm)

        for i in eachindex(test)
            if test[i] < 0.0
                count_dmg = count_dmg + 1;
            elseif test[i] == 0.0
                count_edge = count_edge + 1;
            elseif test[i] > 0.0
                count_norm = count_norm + 1;
                test[i] = 1; #binarize map for presentation
            end
        end
    return test, count_dmg, count_edge, count_norm
end


function perc_dam(test, count_dmg, count_edge, count_norm)

       total_grid = (size(test)[1] * size(test)[2]) - cnt_edge;
       damage_perc = (cnt_dmg / total_grid) * 100

       return damage_perc
end


#Sarcomere Length (μm)
# this is a change
function main_function(filename)

    img = ft_filter(filename)
    rot_img, theta = center_img(img)
    sl = get_sarclength(rot_img)
    sl2 = recenter_img(sl,theta)
    replace!(sl2, NaN=> 0)
    sl2 = sl2[1:1024,1:1024]
    dmg = damage_kernel(sl2)
    count_dmg, count_norm, count_edge = 0,0,0;
    binary_map, cnt_dmg, cnt_edge, cnt_norm = count_dam(dmg, count_dmg, count_edge, count_norm)
    heatmap(binary_map) #binary map (healthy vs. damaged) for presentation
    damage_perc = perc_dam(dmg, cnt_dmg, cnt_edge, cnt_norm)

end

    #save out file to read into MATLAB if necessary
#    writedlm("45%_15OCT19_Series031_location2_1xzoom_ch01.csv", img_processed, ',')


function meanpos(a)
        posMean = mean(a[a .> 0])
        return posMean
end

function remove_shade!(img) #current fxn
    siz = size(img)
    shade = zeros(siz[1], siz[2])
    region = ones(9, 9)./81
    dimg = imfilter(img,centered(region))
    typical = quantile(img[:], 0.3)
    for i in eachindex(img)
        if dimg[i] < typical
        shade[i] = img[i]
        #img[i] *= 1.2
        end
    end
    avg = (meanpos(shade))
    for i in eachindex(shade)
        if shade[i] > 1.2*avg
            img[i] *= 1.5
        end
        if shade[i] < 0.6*avg && shade[i] > 0
            img[i] *= 0.5
        end
    end,
    return img
end
