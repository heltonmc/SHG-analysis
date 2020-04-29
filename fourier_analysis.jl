using Images


#Read image file
function readimg(filename::String)

    #Read file
    img::Array{Float64} = load(filename)
    img = convert(Array{Float64},img)
    img = img./maximum(img)
end

#define window function for fourier transform
function windowfunc(img::Array)

    #=hanwindow = hanning(size(img,1),padding=0,zerophase=false)
    hanwin = hanwindow*hanwindow'
    heatmap(hanwin)=#
    x::Array{Int64,1} = Vector(1:1:size(img,1))
    a::Float64 = size(img,1)/2
    A::Array{Float64,1} = -0.5.*(1 .+ cos.(π.*x./a)) .+ 1

    window::Array{Float64,2} = A*A'
end

#apply padding/window to image and take ft
function FT(image::Array)

    img = image.*windowfunc(image)

    ft = fft(img) |> fftshift
end


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
