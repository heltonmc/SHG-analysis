#here is the OMI analysis file

using FFTW, StatsBase, Plots, ImageView, Images, OffsetArrays, CoordinateTransformations, Dierckx

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

    θ =  mode(round.(angles,digits =0))
    tfm = recenter(RotMatrix(deg2rad(θ)),[size(img)[1]/2,size(img)[1]/2])
    rot_img = warp(img, tfm)
    return rot_img, θ
end


#get sarcomere length estimates from local maximums of columns
function get_sarclength(rot_img)

    sl_data = rot_img
    for in_col in axes(rot_img)[1]

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





#    heatmap(abs.(img_processed))

    #save out file to read into MATLAB if necessary
#    writedlm("45%_15OCT19_Series031_location2_1xzoom_ch01.csv", img_processed, ',')
