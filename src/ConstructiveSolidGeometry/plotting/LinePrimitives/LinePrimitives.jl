include("Line.jl")
include("LineSegment.jl")
include("Arc.jl")

@recipe function f(lines::Array{AbstractLinePrimitive}; n = 30)
    seriescolor --> :orange
    linewidth --> 2
    for line in lines
        @series begin
            label := ""
            line
        end
    end
end
