using ArgParse
using Distributions

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--replicates", "-r"
        help = "Number of replicates per sample"
        arg_type = Int
        default = 2
        "--confidence", "-c"
        help = "Confidence level for intervals"
        arg_type = Int
        default = 95
    end
    return parse_args(s)
end

function main()
    args = parse_commandline()

    (vafs, sampleIDs) = readdlm(STDIN, '\t', header=true)

    @assert size(vafs, 2) % args["replicates"] == 0
    betaDistributions = mapslices(x -> fit(Beta, x), vafs, 2)
    betaDistributions = squeeze(betaDistributions, 2)
    print(size(betaDistributions))
    for betaDist in betaDistributions
        println(betaDist)
    end
    lowerConf = (100 - args["confidence"]) * 0.01 
    lowerCIs = map(x -> quantile(x, lowerConf), betaDistributions)
    println("Lower CIs:")
    for lowerCI in lowerCIs
        println(lowerCI)
    end
    upperConf = args["confidence"] * 0.01 
    upperCIs = map(x -> quantile(x, upperConf), betaDistributions)
    println("Upper CIs:")
    for upperCI in upperCIs
        println(upperCI)
    end
    println("Means:")
    for betaDist in betaDistributions
        println(mean(betaDist))
    end
end

main()
