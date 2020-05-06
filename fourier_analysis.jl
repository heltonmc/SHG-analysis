
"""
    FT_SHG

Standard Library module to load and take fourier transform of Second Harmonic Images
"""

module FT_SHG

using Images,FFTW, DSP, Statistics, ImageFiltering, ImageQualityIndexes

export readimg, windowfunc, ft_image

##### load image #####

"""
    readimg(filename:string)

Load the image file and convert it to a normalized 2D array of Float64 values

# Examples
'''jldoctest
julia> using FT_SHG

julia> readimg("/home/heltonmc/Desktop/Images/Control/11OCT19_Series015_location1_2xzoom_ch01.tif")
1024×1024 Array{Float64,2}:
 0.345455  0.260606  0.224242  0.490909  0.290909 ...
 0.339394  0.309091  0.484848  0.315152  0.424242
 0.181818  0.29697   0.29697   0.345455  0.230303
 0.393939  0.333333  0.224242  0.369697  0.563636
 0.260606  0.424242  0.181818  0.333333  0.248485
 0.193939  0.157576  0.181818  0.309091  0.339394 ...
 ⋮                                             ⋮
'''
"""
function readimg(filename::String)
    img::Array{Float64} = load(filename) .|> Float64
    img = img./maximum(img)
end

##### apply window function to image #####

"""
    windowfunc(img::Array)

Create a 2D window function to apply on an image.

!!! note
    This uses a hanning window explicitly defined. Opportunity to explore other
    window functions. Need to pad image array as it deprecates image
    around borders. Why do we do this?

# Examples
'''jldoctest
julia> using FT_SHG

julia> windowfunc(size(img,1)) # img is your 2D array loaded with readimg function
windowfunc(1024)
1024×1024 Array{Float64,2}:
 8.85925e-11  3.54367e-10  7.97312e-10  1.41741e-9  2.21465e-9 ...
 3.54367e-10  1.41745e-9   3.18922e-9   5.6696e-9   8.8585e-9
 7.97312e-10  3.18922e-9   7.17563e-9   1.27564e-8  1.99313e-8
 1.41741e-9   5.6696e-9    1.27564e-8   2.26775e-8  3.54327e-8
 2.21465e-9   8.8585e-9    1.99313e-8   3.54327e-8  5.5362e-8
 3.18898e-9   1.27558e-8   2.87001e-8   5.10213e-8  7.97185e-8 ...
 ⋮                                                          ⋮
 '''
"""

function windowfunc(pixel_length::Int64)
    #=hanwindow = hanning(size(img,1),padding=0,zerophase=false)
    hanwin = hanwindow*hanwindow'
    heatmap(hanwin)=#
    x::Array{Int64,1} = Vector(1:1:pixel_length)
    a::Float64 = pixel_length/2
    A::Array{Float64,1} = -0.5.*(1 .+ cos.(π.*x./a)) .+ 1

    window::Array{Float64,2} = A*A'
end


##### take fourier transfrom of image #####

"""
    ft_image(image::Array)

Take fourier transfrom of 2D image.

!!! note


# Examples
'''jldoctest
julia> using FT_SHG

julia> ft_image("/home/heltonmc/Desktop/Images/Control/11OCT19_Series015_location1_2xzoom_ch01.tif")
1024×1024 Array{Complex{Float64},2}:
  13.3757+0.0im       4.49254-23.9886im   -27.7867+29.1356im   10.8157-1.39225im ...
  17.4757-26.147im   -9.59154+22.5935im     6.0801-2.1752im    7.16925-6.8986im
 -42.5586+7.61455im   23.0202-11.5096im    15.2907-20.015im   -30.3467+20.3266im
  32.6522+13.1033im  -12.5766+2.97979im   -41.2895+16.8798im   52.9914-12.6014im
 -16.0971-3.28526im   5.28427+9.48263im    26.5788-18.8307im  -26.4525-2.54751im
  5.90529-22.3272im  -17.7353+8.22078im    13.9694+5.24972im  -18.1758+11.5949im ...
         ⋮
 '''
"""
function ft_image(filename::String)
    img::Array{Float64,2} = readimg(filename)
    img = img.*windowfunc(size(img,1)) .|> Float64

    ft = fft(img) |> fftshift
end

#=

function denoise(image)

end
F2 = abs.(I)
noise = maximum(F2[1:100]) # make better noise calculator

for n = 1:size(F2,1)
    for c = 1:size(F2,2)
        if F2[n,c] < 3*noise
            F2[n,c] = 0
        end
    end
end
sum(F2)
end
=#
end


function denoiseimg(img,kernel)

    #kernel = ones(5, 5)./9 # mean filter
    denoised_img = imfilter(img, centered(kernel))
    assess(PSNR(),denoised_img,img)
    assess(SSIM(),denoised_img,img)

    for i in eachindex(denoised_img)
        tol = 0.54
        if abs(denoised_img[i]) < tol
            denoised_img[i] = 0
        else
            denoised_img[i] = 1
        end

    end
    return denoised_img
end



kernel_x = [1 0 -1; 2 0 -2;1 0 -1]
kernel_y = [1 2 1; 0 0 0; -1 -2 -1]
sobel_x = imfilter(img, centered(kernel_x))
grad = imfilter(sobel_x, kernel')
