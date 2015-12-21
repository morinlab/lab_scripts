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
    vafs = readdlm(STDIN, '\t', header=true)

    betaDistributions = fit(Beta, vafs)

    lowerCI = quantile(betaDist, (100 - args.confidence) * 0.1)
    upperCI = quantile(betaDist, (args.confidence * 0.1))
end
