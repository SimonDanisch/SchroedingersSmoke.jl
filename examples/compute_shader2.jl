using GLWindow, GLAbstraction, ModernGL, Reactive, GLFW, GeometryTypes

const window = create_glcontext("Compute Shader", resolution=(512, 512), major=4, minor=5)

# In order to write to a texture, we have to introduce it as image2D.
# local_size_x/y/z layout variables define the work group size.
# gl_GlobalInvocationID is a uvec3 variable giving the global ID of the thread,
# gl_LocalInvocationID is the local index within the work group, and
# gl_WorkGroupID is the work group's index
const shader = comp"""
{{GLSL_VERSION}}

uniform sampler3D velocity;
writeonly uniform image3D result;

uniform vec3 d;

layout (local_size_x = 16, local_size_y = 16, local_size_z = 16) in;

void main() {
    ivec3 res = textureSize(velocity, 0);
    ivec3 xyz = ivec3(gl_GlobalInvocationID.xyz);
    ivec3 im  = ivec3(mod(vec3(xyz-1), vec3(res)));

    float _x = texelFetch(velocity, ivec3(im.x, xyz.yz), 0).x;
    float _y = texelFetch(velocity, ivec3(xyz.x, im.y, xyz.z), 0).y;
    float _z = texelFetch(velocity, ivec3(xyz.xy, im.z), 0).z;

    vec3 v = texelFetch(velocity, xyz, 0).xyz;

    vec3 ff = (v - vec3(_x, _y, _z)) * d;

    imageStore(result, xyz-1, vec4(ff.x+ff.y+ff.z));
}
"""
dims = (256, 256, 256)
const velocity = Texture(Vec3f0, dims)
const f = Texture(Float32, dims)
glBindImageTexture(0, f.id, 0, GL_FALSE, 0, GL_WRITE_ONLY, f.internalformat);

prg = LazyShader(shader)
d = Vec3f0(4)./Vec3f0(dims)
d = d.*d
data = Dict(
    :d    => 1f0/d,
    :velocity => velocity,
)

const robk = RenderObject(data, prg, ()->nothing, ()-> glDispatchCompute(div(256,16), div(256,16), div(256,16)))

render(robk)
glFinish()
gpu_data(f)
