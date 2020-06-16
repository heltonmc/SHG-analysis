#here is the OMI analysis file



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
end


#get index values of local maximums of a column
locmax = findlocalmaxima(img[:,col])
idx = []
for i in eachindex(locmax)
    push!(idx,getindex(locmax,i)[1])
end




    heatmap(abs.(img_processed))

    #save out file to read into MATLAB if necessary
    writedlm("45%_15OCT19_Series031_location2_1xzoom_ch01.csv", img_processed, ',')


#=
