import OpenCL
const cl = OpenCL

import CLFFT
const clfft = CLFFT

include("cl_fft.jl")
readshader(name) = readstring(Pkg.dir("SchroedingersSmoke", "src", name))
const device = first(cl.devices(:cpu))
const ctx    = cl.Context(device)
const queue  = cl.CmdQueue(ctx)
const cl_helpers = readshader("complex.cl") * readshader("helper.cl") 

function compileprogram(ctx, source)
	source = cl_helpers*source
    program = cl.Program(ctx, source=source)
    cl.build!(program, raise=false)
    for (dev, status) in cl.info(program, :build_status)
        dict = cl.info(program, :build_log)
		println(length(split(cl_helpers, "\n")))
        for (k,v) in dict
            println(k)
            println(v)
        end
    end
    program
end


const program_ISF = compileprogram(ctx, readshader("ISF.cl"))
const program_torusDEC = compileprogram(ctx, readshader("torusDEC.cl"))
const program_particle = compileprogram(ctx, readshader("particle.cl"))

macro cl_kernel(name, program, args...)
    kernel = gensym(name)
    quote
        const $kernel = cl.Kernel($program, $(string(name)))
        function $name($(args...))
            cl.call($queue, $kernel, size($(args[1])), nothing, $(args...))
        end
    end
end

const cl_mul = cl.Kernel(program_torusDEC, "cl_mul")
function mul!(a, b)
	cl.call(queue, cl_mul, length(a), nothing, a, b)
end


const cl_mul_cf = cl.Kernel(program_torusDEC, "cl_mul_cf")
function mul!(a::cl.CLArray{Complex}, b::cl.CLArray{Complex})
	cl.call(queue, cl_mul_cf, length(a), nothing, a, b)
end