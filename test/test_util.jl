# https://discourse.julialang.org/t/skipping-a-whole-testset/65006/4
using Test: DefaultTestSet, Broken, parse_testset_args, Test, finish
macro testset_skip(args...)
    isempty(args) && error("No arguments to @testset_skip")
    length(args) < 2 && error("First argument to @testset_skip giving reason for "
          *
          "skipping is required")
    skip_reason = args[1]
    desc, testsettype, options = parse_testset_args(args[2:(end - 1)])
    ex = quote
        # record the reason for the skip in the description, and mark the tests as
        # broken, but don't run tests
        local ts = DefaultTestSet(string($desc, " - ", $skip_reason))
        push!(ts.results, Broken(:skipped, "skipped tests"))
        local ret = finish(ts)
        ret
    end
    return ex
end
