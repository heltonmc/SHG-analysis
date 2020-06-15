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
    heatmap(abs.(img_processed))

    #save out file to read into MATLAB if necessary
    writedlm("45%_15OCT19_Series031_location2_1xzoom_ch01.csv", img_processed, ',')


#=
