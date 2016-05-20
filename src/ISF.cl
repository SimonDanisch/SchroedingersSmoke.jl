typedef struct _ISF{
  float hbar;
  float dt;
}ISF;

typedef struct _Torus{
  float hbar;
  float dt;
}Torus;



"""
extracts velocity 1-form from (psi1,psi2).
If hbar argument is empty, hbar=1 is assumed.
"""
void VelocityOneForm(
        Torus obj,
        __global const float2* psi1,
        __global const float2* psi2,
        __global float3* velocity,
        float hbar
    ){
    int x = get_global_id(0);
    int y = get_global_id(1);
    int z = get_global_id(2);
    ixp = mod(x, obj.resx) + 1;
    iyp = mod(y, obj.resy) + 1;
    izp = mod(z, obj.resz) + 1;
    psi1c = conj(psi1[x,y,z]); psi2c = conj(psi1[x,y,z]);
    vx = angle(psi1c* psi1[ixp,y,z] + psi2c* psi2[ixp,y,z]);
    vy = angle(psi1c.*psi1[x,iyp,z] + psi2c.*psi2[x,iyp,z]);
    vz = angle(psi1c.*psi1[x,y,izp] + psi2c.*psi2[x,y,izp]);
    velocity[x,y,z] = float3(vx*hbar, vy*hbar, vz*hbar);
}


"""
multiplies exp(i*q) to (psi1,psi2)
"""
void GaugeTransform(
        __global const float2* psi1,
        __global const float2* psi2,
        __global float2* q,
    )
    int i = get_global_id(0);
    float eiq = exp(vec2(0,1)*q[i]);
    psi1[i] *= eiq;
    psi2[i] *= eiq;
end



void Hopf(
        __global const float2* psi1,
        __global const float2* psi2,
        __global float3* s,
    ){
    int i = get_global_id(0);
    s[i] = _hopf(psi1[i], psi2[i]);
}

float3 _hopf(float2 psi1, float2 psi2)
    a = real(psi1);
    b = imag(psi1);
    c = real(psi2);
    d = imag(psi2);
    sx = 2*(a*c + b*d);
    sy = 2*(a*d - b*c);
    sz = a^2 + b^2 - c^2 - d^2;
    float3(sx,sy,sz);
end

float2 _normalize(float2 psi1, float2 psi2){
    float2(psi1/norm, psi2/norm);
}

void Normalize(
        __global float2* psi1,
        __global float2* psi2
    ){
    int i = get_global_id(0);
    float psi1_ = psi1[i];
    float psi2_ = psi1[i];
    float norm = sqrt(abs(psi1_)^2 + abs(psi2_)^2);
    psi1[i] = psi1_/norm;
    psi2[i] = psi2_/norm;
}
