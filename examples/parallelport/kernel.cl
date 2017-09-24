 // dependencies
 // CLArrays.HostPtr{Tuple{Float32,Float32,Float32}}
 struct  __attribute__ ((packed)) TYPHostPtr_float3{
     int ptr;
 };
 typedef struct TYPHostPtr_float3 HostPtr_float3;

 // CLArrays.DeviceArray{Tuple{Float32,Float32,Float32},3,CLArrays.HostPtr{Tuple{Float32,Float32,Float32}}}
 struct  __attribute__ ((packed)) TYPDeviceArray_float3_3_HostPtr_float3{
     HostPtr_float3 ptr;
     uint3 size;
 };
 typedef struct TYPDeviceArray_float3_3_HostPtr_float3 DeviceArray_float3_3_HostPtr_float3;

 // CLArrays.DeviceArray{Tuple{Float32,Float32,Float32},3,Transpiler.CLIntrinsics.GlobalPointer{Tuple{Float32,Float32,Float32}}}
 struct  __attribute__ ((packed)) TYPDeviceArray_float3_3___global1float3121{
     __global float3 *  ptr;
     uint3 size;
 };
 typedef struct TYPDeviceArray_float3_3___global1float3121 DeviceArray_float3_3___global1float3121;

 // CLArrays.#reconstruct
 __constant int FUNC_INST_CLArrays34reconstruct = 0;
 typedef int CLArrays34reconstruct; // empty type emitted as an int
 // Any
 typedef int Any; // placeholder type instance
 __constant Any TYP_INST_Any = 0;

 // Type{Transpiler.CLIntrinsics.GlobalPointer{Tuple{Float32,Float32,Float32}}}
 typedef int Type5Transpiler3CLIntrinsics3GlobalPointer5Tuple5Float327Float327Float32666; // placeholder type instance
 __constant Type5Transpiler3CLIntrinsics3GlobalPointer5Tuple5Float327Float327Float32666 TYP_INST_Type5Transpiler3CLIntrinsics3GlobalPointer5Tuple5Float327Float327Float32666 = 0;

 // (CLArrays.DeviceArray{Tuple{Float32,Float32,Float32},3,Transpiler.CLIntrinsics.GlobalPointer{Tuple{Float32,Float32,Float32}}}, Tuple{Transpiler.CLIntrinsics.GlobalPointer{Tuple{Float32,Float32,Float32}},Tuple{UInt32,UInt32,UInt32}})
 DeviceArray_float3_3___global1float3121 x8DeviceArray_float3_3___global1float31219_33(__global float3 *  ptr, uint3 size)
 {
     return (DeviceArray_float3_3___global1float3121){ptr, size};
 }
 // Symbol

 // (CLArrays.reconstruct, Tuple{CLArrays.DeviceArray{Tuple{Float32,Float32,Float32},3,CLArrays.HostPtr{Tuple{Float32,Float32,Float32}}},Transpiler.CLIntrinsics.GlobalPointer{Tuple{Float32,Float32,Float32}}})
 DeviceArray_float3_3___global1float3121 reconstruct_34(DeviceArray_float3_3_HostPtr_float3 x, __global float3 *  ptr)
 {
     return x8DeviceArray_float3_3___global1float31219_33(ptr, x.size);
 }
 // #staggered_velocity
 __constant int FUNC_INST_x4staggered_velocity = 0;
 typedef int x4staggered_velocity; // empty type emitted as an int
 // ##19#23
 __constant int FUNC_INST_x4419423 = 0;
 typedef int x4419423; // empty type emitted as an int
 // ##20#24
 __constant int FUNC_INST_x4420424 = 0;
 typedef int x4420424; // empty type emitted as an int
 // ##21#25
 __constant int FUNC_INST_x4421425 = 0;
 typedef int x4421425; // empty type emitted as an int
 // ##22#26
 __constant int FUNC_INST_x4422426 = 0;
 typedef int x4422426; // empty type emitted as an int
 // Base.#mod
 __constant int FUNC_INST_Base34mod = 0;
 typedef int Base34mod; // empty type emitted as an int
 // Base.#broadcast
 __constant int FUNC_INST_Base34broadcast = 0;
 typedef int Base34broadcast; // empty type emitted as an int
 // Tuple{Tuple{UInt32,UInt32,UInt32}}
 struct  __attribute__ ((packed)) TYPTuple_uint3{
     uint3 field1;
 };
 typedef struct TYPTuple_uint3 Tuple_uint3;

 // Base.#map
 __constant int FUNC_INST_Base34map = 0;
 typedef int Base34map; // empty type emitted as an int
 // Base.#promote
 __constant int FUNC_INST_Base34promote = 0;
 typedef int Base34promote; // empty type emitted as an int
 // Type{Float32}
 typedef int Type5Float326; // placeholder type instance
 __constant Type5Float326 TYP_INST_Type5Float326 = 0;

 // Type{UInt32}
 typedef int Type5UInt326; // placeholder type instance
 __constant Type5UInt326 TYP_INST_Type5UInt326 = 0;

 // (promote, Tuple{Float32,UInt32})
 float2 promote_154(float x, uint y)
 {
     return (float2){x, (float){y}};
 }
 // (mod, Tuple{Float32,Float32})
 float mod_24(float x, float y)
 {
     float r;
     r = remainder(x, y);
     if(r == 0){
         return copysign(r, y);
     };
     if((r > 0) ^ (y > 0)){
         return r + y;
     };
     return r;
 }
 // (mod, Tuple{Float32,UInt32})
 float mod_154(float x, uint y)
 {
     float2 x44_apply_tmp41504;
     x44_apply_tmp41504 = promote_154(x, y);
     return mod_24(x44_apply_tmp41504.s0, x44_apply_tmp41504.s1);
 }
 // (map, Tuple{Base.#mod,Tuple{Float32,Float32},Tuple{UInt32,UInt32}})
 float2 map_155(Base34mod f, float2 t, uint2 s)
 {
     return (float2){mod_154(t.s0, s.s0), mod_154(t.s1, s.s1)};
 }
 // Base.#tail
 __constant int FUNC_INST_Base34tail = 0;
 typedef int Base34tail; // empty type emitted as an int
 // Base.#argtail
 __constant int FUNC_INST_Base34argtail = 0;
 typedef int Base34argtail; // empty type emitted as an int
 // (Base.argtail, Tuple{Float32,Float32,Float32})
 float2 argtail_156(float x, float2 rest)
 {
     return rest;
 }
 // (Base.tail, Tuple{Tuple{Float32,Float32,Float32}})
 float2 tail_122(float3 x)
 {
     float3 x44_apply_tmp41505;
     x44_apply_tmp41505 = x;
     return argtail_156(x44_apply_tmp41505.s0, (float2){x44_apply_tmp41505.s1, x44_apply_tmp41505.s2});
 }
 // (Base.argtail, Tuple{UInt32,UInt32,UInt32})
 uint2 argtail_8(uint x, uint2 rest)
 {
     return rest;
 }
 // (Base.tail, Tuple{Tuple{UInt32,UInt32,UInt32}})
 uint2 tail_9(uint3 x)
 {
     uint3 x44_apply_tmp41506;
     x44_apply_tmp41506 = x;
     return argtail_8(x44_apply_tmp41506.s0, (uint2){x44_apply_tmp41506.s1, x44_apply_tmp41506.s2});
 }
 // (map, Tuple{Base.#mod,Tuple{Float32,Float32,Float32},Tuple{UInt32,UInt32,UInt32}})
 float3 map_157(Base34mod f, float3 t, uint3 s)
 {
     float2 x44_apply_tmp41503;
     x44_apply_tmp41503 = map_155(f, tail_122(t), tail_9(s));
     return (float3){mod_154(t.s0, s.s0), x44_apply_tmp41503.s0, x44_apply_tmp41503.s1};
 }
 // (broadcast, Tuple{Base.#mod,Tuple{Float32,Float32,Float32},Tuple{UInt32,UInt32,UInt32}})
 float3 broadcast_157(Base34mod f, float3 t, Tuple_uint3 ts)
 {
     Tuple_uint3 x44_apply_tmp41502;
     x44_apply_tmp41502 = ts;
     return map_157(f, t, x44_apply_tmp41502.field1);
 }
 // Type{##19#23}
 typedef int Type544194236; // placeholder type instance
 __constant Type544194236 TYP_INST_Type544194236 = 0;

 // Tuple{Tuple{Float32,Float32,Float32},UInt32}
 struct  __attribute__ ((packed)) TYPTuple_float3_uint{
     float3 field1;
     uint field2;
 };
 typedef struct TYPTuple_float3_uint Tuple_float3_uint;

 // Tuple
 typedef int EmptyTuple; // placeholder type instance
 __constant EmptyTuple TYP_INST_EmptyTuple = 0;

 // Type{Tuple}
 typedef int Type5Tuple6; // placeholder type instance
 __constant Type5Tuple6 TYP_INST_Type5Tuple6 = 0;

 // Base.Broadcast.#containertype
 __constant int FUNC_INST_Base3Broadcast34containertype = 0;
 typedef int Base3Broadcast34containertype; // empty type emitted as an int
 // Base.Broadcast.#promote_containertype
 __constant int FUNC_INST_Base3Broadcast34promote_containertype = 0;
 typedef int Base3Broadcast34promote_containertype; // empty type emitted as an int
 // (Base.Broadcast.promote_containertype, Tuple{Type{Tuple},Type{Tuple}})
 Type5Tuple6 promote_containertype_92(Type5Tuple6 x4unused4, Type5Tuple6 x4unused41)
 {
     return TYP_INST_Type5Tuple6;
 }
 // Type{Tuple{Float32,Float32,Float32}}
 typedef int Type5Tuple5Float327Float327Float3266; // placeholder type instance
 __constant Type5Tuple5Float327Float327Float3266 TYP_INST_Type5Tuple5Float327Float327Float3266 = 0;

 // Base.Broadcast.#_containertype
 __constant int FUNC_INST_Base3Broadcast34_containertype = 0;
 typedef int Base3Broadcast34_containertype; // empty type emitted as an int
 // (Base.Broadcast._containertype, Tuple{Type{Tuple{Float32,Float32,Float32}}})
 Type5Tuple6 _containertype_158(Type5Tuple5Float327Float327Float3266 x4unused4)
 {
     return TYP_INST_Type5Tuple6;
 }
 // (Base.Broadcast.containertype, Tuple{Tuple{Float32,Float32,Float32}})
 Type5Tuple6 containertype_122(float3 x)
 {
     return _containertype_158(TYP_INST_Type5Tuple5Float327Float327Float3266);
 }
 // Type{Any}
 typedef int Type5Any6; // placeholder type instance
 __constant Type5Any6 TYP_INST_Type5Any6 = 0;

 // (Base.Broadcast.promote_containertype, Tuple{Type{Tuple},Type{Any}})
 Type5Tuple6 promote_containertype_55(Type5Tuple6 x4unused4, Type5Any6 x4unused41)
 {
     return TYP_INST_Type5Tuple6;
 }
 // (Base.Broadcast._containertype, Tuple{Type{UInt32}})
 Type5Any6 _containertype_159(Type5UInt326 x4unused4)
 {
     return TYP_INST_Type5Any6;
 }
 // (Base.Broadcast.containertype, Tuple{UInt32})
 Type5Any6 containertype_11(uint x)
 {
     return _containertype_159(TYP_INST_Type5UInt326);
 }
 // (Base.Broadcast.containertype, Tuple{Tuple{Float32,Float32,Float32},UInt32})
 Type5Tuple6 containertype_160(float3 ct1, uint ct2)
 {
     return promote_containertype_55(containertype_122(ct1), containertype_11(ct2));
 }
 // (Base.Broadcast.containertype, Tuple{Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},UInt32})
 Type5Tuple6 containertype_161(float3 ct1, float3 ct2, uint cts)
 {
     uint x44_apply_tmp41509;
     x44_apply_tmp41509 = cts;
     return promote_containertype_92(containertype_122(ct1), containertype_160(ct2, x44_apply_tmp41509));
 }
 // Base.Broadcast.#broadcast_c
 __constant int FUNC_INST_Base3Broadcast34broadcast_c = 0;
 typedef int Base3Broadcast34broadcast_c; // empty type emitted as an int
 // Base.Broadcast.#first_tuple
 __constant int FUNC_INST_Base3Broadcast34first_tuple = 0;
 typedef int Base3Broadcast34first_tuple; // empty type emitted as an int
 // (Base.Broadcast.first_tuple, Tuple{Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},UInt32})
 float3 first_tuple_161(float3 A, Tuple_float3_uint Bs)
 {
     return A;
 }
 // Base.Broadcast.#tuplebroadcast
 __constant int FUNC_INST_Base3Broadcast34tuplebroadcast = 0;
 typedef int Base3Broadcast34tuplebroadcast; // empty type emitted as an int
 // Tuple{Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},UInt32}
 struct  __attribute__ ((packed)) TYPTuple_float3_float3_uint{
     float3 field1;
     float3 field2;
     uint field3;
 };
 typedef struct TYPTuple_float3_float3_uint Tuple_float3_float3_uint;

 // Base.Broadcast.##3#4{##19#23,Tuple{Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},UInt32}}
 struct  __attribute__ ((packed)) TYPBase3Broadcast344344544194237Tuple5Tuple5Float327Float327Float3267Tuple5Float327Float327Float3267UInt3266{
     x4419423 f;
     Tuple_float3_float3_uint As;
 };
 typedef struct TYPBase3Broadcast344344544194237Tuple5Tuple5Float327Float327Float3267Tuple5Float327Float327Float3267UInt3266 Base3Broadcast344344544194237Tuple5Tuple5Float327Float327Float3267Tuple5Float327Float327Float3267UInt3266;

 // Val{3}
 typedef int Val_3; // empty type emitted as an int
 // Type{Val{3}}
 typedef int Type5Val5366; // placeholder type instance
 __constant Type5Val5366 TYP_INST_Type5Val5366 = 0;

 // Base.#ntuple
 __constant int FUNC_INST_Base34ntuple = 0;
 typedef int Base34ntuple; // empty type emitted as an int
 // Tuple{Float32,Float32,UInt32}
 struct  __attribute__ ((packed)) TYPTuple_float_float_uint{
     float field1;
     float field2;
     uint field3;
 };
 typedef struct TYPTuple_float_float_uint Tuple_float_float_uint;

 // Base.Broadcast.#tuplebroadcast_getargs
 __constant int FUNC_INST_Base3Broadcast34tuplebroadcast_getargs = 0;
 typedef int Base3Broadcast34tuplebroadcast_getargs; // empty type emitted as an int
 // Base.Broadcast.#_broadcast_getindex
 __constant int FUNC_INST_Base3Broadcast34_broadcast_getindex = 0;
 typedef int Base3Broadcast34_broadcast_getindex; // empty type emitted as an int
 // (Base.Broadcast._broadcast_getindex, Tuple{Type{Tuple},Tuple{Float32,Float32,Float32},Int64})
 float _broadcast_getindex_162(Any x4unused4, float3 A, long I)
 {
     return ((float*)&A)[I - 1];
 }
 // (Base.Broadcast._broadcast_getindex, Tuple{Tuple{Float32,Float32,Float32},Int64})
 float _broadcast_getindex_163(float3 A, long I)
 {
     return _broadcast_getindex_162(containertype_122(A), A, I);
 }
 // Base.#first
 __constant int FUNC_INST_Base34first = 0;
 typedef int Base34first; // empty type emitted as an int
 // (first, Tuple{Tuple{Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},UInt32}})
 float3 first_164(Tuple_float3_float3_uint t)
 {
     return t.field1;
 }
 // Tuple{Float32,UInt32}
 struct  __attribute__ ((packed)) TYPTuple_float_uint{
     float field1;
     uint field2;
 };
 typedef struct TYPTuple_float_uint Tuple_float_uint;

 // (first, Tuple{Tuple{Tuple{Float32,Float32,Float32},UInt32}})
 float3 first_165(Tuple_float3_uint t)
 {
     return t.field1;
 }
 // (Base.Broadcast._broadcast_getindex, Tuple{Type{Any},UInt32,Int64})
 uint _broadcast_getindex_166(Any x4unused4, uint A, long I)
 {
     return A;
 }
 // (Base.Broadcast._broadcast_getindex, Tuple{UInt32,Int64})
 uint _broadcast_getindex_167(uint A, long I)
 {
     return _broadcast_getindex_166(containertype_11(A), A, I);
 }
 // (first, Tuple{Tuple{UInt32}})
 uint first_12(uint t)
 {
     return t;
 }
 // Tuple{}
 typedef int EmptyTuple_; // empty type emitted as an int
 // (Base.Broadcast.tuplebroadcast_getargs, Tuple{Tuple{},Int64})
 EmptyTuple_ tuplebroadcast_getargs_66(EmptyTuple_ x4unused4, long k)
 {
     return (EmptyTuple_){0.0f};
 }
 // (Base.argtail, Tuple{UInt32})
 EmptyTuple_ argtail_11(uint x, EmptyTuple_ rest)
 {
     return rest;
 }
 // (Base.tail, Tuple{Tuple{UInt32}})
 EmptyTuple_ tail_12(uint x)
 {
     uint x44_apply_tmp41516;
     x44_apply_tmp41516 = x;
     return argtail_11(x44_apply_tmp41516, (EmptyTuple_){0.0f});
 }
 // (Base.Broadcast.tuplebroadcast_getargs, Tuple{Tuple{UInt32},Int64})
 uint tuplebroadcast_getargs_168(uint As, long k)
 {
     EmptyTuple_ x44_apply_tmp41515;
     x44_apply_tmp41515 = tuplebroadcast_getargs_66(tail_12(As), k);
     return (uint){_broadcast_getindex_167(first_12(As), k)};
 }
 // (Base.argtail, Tuple{Tuple{Float32,Float32,Float32},UInt32})
 uint argtail_160(float3 x, uint rest)
 {
     return rest;
 }
 // (Base.tail, Tuple{Tuple{Tuple{Float32,Float32,Float32},UInt32}})
 uint tail_165(Tuple_float3_uint x)
 {
     Tuple_float3_uint x44_apply_tmp41517;
     x44_apply_tmp41517 = x;
     return argtail_160(x44_apply_tmp41517.field1, (uint){x44_apply_tmp41517.field2});
 }
 // (Base.Broadcast.tuplebroadcast_getargs, Tuple{Tuple{Tuple{Float32,Float32,Float32},UInt32},Int64})
 Tuple_float_uint tuplebroadcast_getargs_169(Tuple_float3_uint As, long k)
 {
     uint x44_apply_tmp41514;
     x44_apply_tmp41514 = tuplebroadcast_getargs_168(tail_165(As), k);
     return (Tuple_float_uint){_broadcast_getindex_163(first_165(As), k), x44_apply_tmp41514};
 }
 // (Base.argtail, Tuple{Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},UInt32})
 Tuple_float3_uint argtail_161(float3 x, Tuple_float3_uint rest)
 {
     return rest;
 }
 // (Base.tail, Tuple{Tuple{Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},UInt32}})
 Tuple_float3_uint tail_164(Tuple_float3_float3_uint x)
 {
     Tuple_float3_float3_uint x44_apply_tmp41518;
     x44_apply_tmp41518 = x;
     return argtail_161(x44_apply_tmp41518.field1, (Tuple_float3_uint){x44_apply_tmp41518.field2, x44_apply_tmp41518.field3});
 }
 // (Base.Broadcast.tuplebroadcast_getargs, Tuple{Tuple{Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},UInt32},Int64})
 Tuple_float_float_uint tuplebroadcast_getargs_170(Tuple_float3_float3_uint As, long k)
 {
     Tuple_float_uint x44_apply_tmp41513;
     x44_apply_tmp41513 = tuplebroadcast_getargs_169(tail_164(As), k);
     return (Tuple_float_float_uint){_broadcast_getindex_163(first_164(As), k), x44_apply_tmp41513.field1, x44_apply_tmp41513.field2};
 }
 // (#19, Tuple{Float32,Float32,UInt32})
 uint x419_171(float x4temp4, float x4temp41, uint x4temp42)
 {
     uint x4temp45;
     float x4temp44;
     float x4temp43;
     x4temp43 = x4temp4 / x4temp41;
     x4temp44 = floor(x4temp43);
     x4temp45 = (uint)(x4temp44);
     return x4temp45 + x4temp42;
 }
 // (Base.Broadcast.##3#4{##19#23,Tuple{Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},UInt32}}, Tuple{Int64})
 uint x8Base3Broadcast344344544194237Tuple5Tuple5Float327Float327Float3267Tuple5Float327Float327Float3267UInt32669_70(Base3Broadcast344344544194237Tuple5Tuple5Float327Float327Float3267Tuple5Float327Float327Float3267UInt3266 x4self4, long k)
 {
     Tuple_float_float_uint x44_apply_tmp41512;
     x44_apply_tmp41512 = tuplebroadcast_getargs_170(x4self4.As, k);
     return x419_171(x44_apply_tmp41512.field1, x44_apply_tmp41512.field2, x44_apply_tmp41512.field3);
 }
 // (ntuple, Tuple{Base.Broadcast.##3#4{##19#23,Tuple{Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},UInt32}},Type{Val{3}}})
 uint3 ntuple_172(Base3Broadcast344344544194237Tuple5Tuple5Float327Float327Float3267Tuple5Float327Float327Float3267UInt3266 f, Type5Val5366 x4unused4)
 {
     return (uint3){x8Base3Broadcast344344544194237Tuple5Tuple5Float327Float327Float3267Tuple5Float327Float327Float3267UInt32669_70(f, 1), x8Base3Broadcast344344544194237Tuple5Tuple5Float327Float327Float3267Tuple5Float327Float327Float3267UInt32669_70(f, 2), x8Base3Broadcast344344544194237Tuple5Tuple5Float327Float327Float3267Tuple5Float327Float327Float3267UInt32669_70(f, 3)};
 }
 // (Base.Broadcast.tuplebroadcast, Tuple{##19#23,Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},UInt32})
 uint3 tuplebroadcast_173(x4419423 f, float3 x4unused4, Tuple_float3_float3_uint As)
 {
     Base3Broadcast344344544194237Tuple5Tuple5Float327Float327Float3267Tuple5Float327Float327Float3267UInt3266 x43;
     x43 = (Base3Broadcast344344544194237Tuple5Tuple5Float327Float327Float3267Tuple5Float327Float327Float3267UInt3266){f, As};
     Base3Broadcast344344544194237Tuple5Tuple5Float327Float327Float3267Tuple5Float327Float327Float3267UInt3266 _ssavalue_0;
     _ssavalue_0 = x43;
     Type5Val5366 _ssavalue_1;
     _ssavalue_1 = TYP_INST_Type5Val5366;
     return ntuple_172(_ssavalue_0, _ssavalue_1);
 }
 // (Base.Broadcast.broadcast_c, Tuple{##19#23,Type{Tuple},Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},UInt32})
 uint3 broadcast_c_174(x4419423 f, Any x4unused4, float3 A, Tuple_float3_uint Bs)
 {
     Tuple_float3_uint x44_apply_tmp41510;
     x44_apply_tmp41510 = Bs;
     Tuple_float3_uint x44_apply_tmp41511;
     x44_apply_tmp41511 = Bs;
     return tuplebroadcast_173(f, first_tuple_161(A, (Tuple_float3_uint){x44_apply_tmp41510.field1, x44_apply_tmp41510.field2}), (Tuple_float3_float3_uint){A, x44_apply_tmp41511.field1, x44_apply_tmp41511.field2});
 }
 // (broadcast, Tuple{##19#23,Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},UInt32})
 uint3 broadcast_175(x4419423 f, float3 A, Tuple_float3_uint Bs)
 {
     Tuple_float3_uint x44_apply_tmp41507;
     x44_apply_tmp41507 = Bs;
     Tuple_float3_uint x44_apply_tmp41508;
     x44_apply_tmp41508 = Bs;
     return broadcast_c_174(f, containertype_161(A, x44_apply_tmp41507.field1, (uint){x44_apply_tmp41507.field2}), A, (Tuple_float3_uint){x44_apply_tmp41508.field1, x44_apply_tmp41508.field2});
 }
 // Type{##20#24}
 typedef int Type544204246; // placeholder type instance
 __constant Type544204246 TYP_INST_Type544204246 = 0;

 // Tuple{Tuple{UInt32,UInt32,UInt32},UInt32}
 struct  __attribute__ ((packed)) TYPTuple_uint3_uint{
     uint3 field1;
     uint field2;
 };
 typedef struct TYPTuple_uint3_uint Tuple_uint3_uint;

 // Type{Tuple{UInt32,UInt32,UInt32}}
 typedef int Type5Tuple5UInt327UInt327UInt3266; // placeholder type instance
 __constant Type5Tuple5UInt327UInt327UInt3266 TYP_INST_Type5Tuple5UInt327UInt327UInt3266 = 0;

 // (Base.Broadcast._containertype, Tuple{Type{Tuple{UInt32,UInt32,UInt32}}})
 Type5Tuple6 _containertype_176(Type5Tuple5UInt327UInt327UInt3266 x4unused4)
 {
     return TYP_INST_Type5Tuple6;
 }
 // (Base.Broadcast.containertype, Tuple{Tuple{UInt32,UInt32,UInt32}})
 Type5Tuple6 containertype_9(uint3 x)
 {
     return _containertype_176(TYP_INST_Type5Tuple5UInt327UInt327UInt3266);
 }
 // (Base.Broadcast.containertype, Tuple{Tuple{UInt32,UInt32,UInt32},UInt32})
 Type5Tuple6 containertype_10(uint3 ct1, uint ct2)
 {
     return promote_containertype_55(containertype_9(ct1), containertype_11(ct2));
 }
 // (Base.Broadcast.containertype, Tuple{Tuple{UInt32,UInt32,UInt32},Tuple{UInt32,UInt32,UInt32},UInt32})
 Type5Tuple6 containertype_177(uint3 ct1, uint3 ct2, uint cts)
 {
     uint x44_apply_tmp41521;
     x44_apply_tmp41521 = cts;
     return promote_containertype_92(containertype_9(ct1), containertype_10(ct2, x44_apply_tmp41521));
 }
 // (Base.Broadcast.first_tuple, Tuple{Tuple{UInt32,UInt32,UInt32},Tuple{UInt32,UInt32,UInt32},UInt32})
 uint3 first_tuple_177(uint3 A, Tuple_uint3_uint Bs)
 {
     return A;
 }
 // Tuple{Tuple{UInt32,UInt32,UInt32},Tuple{UInt32,UInt32,UInt32},UInt32}
 struct  __attribute__ ((packed)) TYPTuple_uint3_uint3_uint{
     uint3 field1;
     uint3 field2;
     uint field3;
 };
 typedef struct TYPTuple_uint3_uint3_uint Tuple_uint3_uint3_uint;

 // Base.Broadcast.##3#4{##20#24,Tuple{Tuple{UInt32,UInt32,UInt32},Tuple{UInt32,UInt32,UInt32},UInt32}}
 struct  __attribute__ ((packed)) TYPBase3Broadcast344344544204247Tuple5Tuple5UInt327UInt327UInt3267Tuple5UInt327UInt327UInt3267UInt3266{
     x4420424 f;
     Tuple_uint3_uint3_uint As;
 };
 typedef struct TYPBase3Broadcast344344544204247Tuple5Tuple5UInt327UInt327UInt3267Tuple5UInt327UInt327UInt3267UInt3266 Base3Broadcast344344544204247Tuple5Tuple5UInt327UInt327UInt3267Tuple5UInt327UInt327UInt3267UInt3266;

 // (Base.Broadcast._broadcast_getindex, Tuple{Type{Tuple},Tuple{UInt32,UInt32,UInt32},Int64})
 uint _broadcast_getindex_178(Any x4unused4, uint3 A, long I)
 {
     return ((uint*)&A)[I - 1];
 }
 // (Base.Broadcast._broadcast_getindex, Tuple{Tuple{UInt32,UInt32,UInt32},Int64})
 uint _broadcast_getindex_179(uint3 A, long I)
 {
     return _broadcast_getindex_178(containertype_9(A), A, I);
 }
 // (first, Tuple{Tuple{Tuple{UInt32,UInt32,UInt32},Tuple{UInt32,UInt32,UInt32},UInt32}})
 uint3 first_180(Tuple_uint3_uint3_uint t)
 {
     return t.field1;
 }
 // (first, Tuple{Tuple{Tuple{UInt32,UInt32,UInt32},UInt32}})
 uint3 first_181(Tuple_uint3_uint t)
 {
     return t.field1;
 }
 // (Base.argtail, Tuple{Tuple{UInt32,UInt32,UInt32},UInt32})
 uint argtail_10(uint3 x, uint rest)
 {
     return rest;
 }
 // (Base.tail, Tuple{Tuple{Tuple{UInt32,UInt32,UInt32},UInt32}})
 uint tail_181(Tuple_uint3_uint x)
 {
     Tuple_uint3_uint x44_apply_tmp41527;
     x44_apply_tmp41527 = x;
     return argtail_10(x44_apply_tmp41527.field1, (uint){x44_apply_tmp41527.field2});
 }
 // (Base.Broadcast.tuplebroadcast_getargs, Tuple{Tuple{Tuple{UInt32,UInt32,UInt32},UInt32},Int64})
 uint2 tuplebroadcast_getargs_182(Tuple_uint3_uint As, long k)
 {
     uint x44_apply_tmp41526;
     x44_apply_tmp41526 = tuplebroadcast_getargs_168(tail_181(As), k);
     return (uint2){_broadcast_getindex_179(first_181(As), k), x44_apply_tmp41526};
 }
 // (Base.argtail, Tuple{Tuple{UInt32,UInt32,UInt32},Tuple{UInt32,UInt32,UInt32},UInt32})
 Tuple_uint3_uint argtail_177(uint3 x, Tuple_uint3_uint rest)
 {
     return rest;
 }
 // (Base.tail, Tuple{Tuple{Tuple{UInt32,UInt32,UInt32},Tuple{UInt32,UInt32,UInt32},UInt32}})
 Tuple_uint3_uint tail_180(Tuple_uint3_uint3_uint x)
 {
     Tuple_uint3_uint3_uint x44_apply_tmp41528;
     x44_apply_tmp41528 = x;
     return argtail_177(x44_apply_tmp41528.field1, (Tuple_uint3_uint){x44_apply_tmp41528.field2, x44_apply_tmp41528.field3});
 }
 // (Base.Broadcast.tuplebroadcast_getargs, Tuple{Tuple{Tuple{UInt32,UInt32,UInt32},Tuple{UInt32,UInt32,UInt32},UInt32},Int64})
 uint3 tuplebroadcast_getargs_183(Tuple_uint3_uint3_uint As, long k)
 {
     uint2 x44_apply_tmp41525;
     x44_apply_tmp41525 = tuplebroadcast_getargs_182(tail_180(As), k);
     return (uint3){_broadcast_getindex_179(first_180(As), k), x44_apply_tmp41525.s0, x44_apply_tmp41525.s1};
 }
 // (#20, Tuple{UInt32,UInt32,UInt32})
 uint x420_8(uint x4temp4, uint x4temp41, uint x4temp42)
 {
     uint x4temp43;
     x4temp43 = x4temp4 % x4temp41;
     return x4temp43 + x4temp42;
 }
 // (Base.Broadcast.##3#4{##20#24,Tuple{Tuple{UInt32,UInt32,UInt32},Tuple{UInt32,UInt32,UInt32},UInt32}}, Tuple{Int64})
 uint x8Base3Broadcast344344544204247Tuple5Tuple5UInt327UInt327UInt3267Tuple5UInt327UInt327UInt3267UInt32669_70(Base3Broadcast344344544204247Tuple5Tuple5UInt327UInt327UInt3267Tuple5UInt327UInt327UInt3267UInt3266 x4self4, long k)
 {
     uint3 x44_apply_tmp41524;
     x44_apply_tmp41524 = tuplebroadcast_getargs_183(x4self4.As, k);
     return x420_8(x44_apply_tmp41524.s0, x44_apply_tmp41524.s1, x44_apply_tmp41524.s2);
 }
 // (ntuple, Tuple{Base.Broadcast.##3#4{##20#24,Tuple{Tuple{UInt32,UInt32,UInt32},Tuple{UInt32,UInt32,UInt32},UInt32}},Type{Val{3}}})
 uint3 ntuple_184(Base3Broadcast344344544204247Tuple5Tuple5UInt327UInt327UInt3267Tuple5UInt327UInt327UInt3267UInt3266 f, Type5Val5366 x4unused4)
 {
     return (uint3){x8Base3Broadcast344344544204247Tuple5Tuple5UInt327UInt327UInt3267Tuple5UInt327UInt327UInt3267UInt32669_70(f, 1), x8Base3Broadcast344344544204247Tuple5Tuple5UInt327UInt327UInt3267Tuple5UInt327UInt327UInt3267UInt32669_70(f, 2), x8Base3Broadcast344344544204247Tuple5Tuple5UInt327UInt327UInt3267Tuple5UInt327UInt327UInt3267UInt32669_70(f, 3)};
 }
 // (Base.Broadcast.tuplebroadcast, Tuple{##20#24,Tuple{UInt32,UInt32,UInt32},Tuple{UInt32,UInt32,UInt32},Tuple{UInt32,UInt32,UInt32},UInt32})
 uint3 tuplebroadcast_185(x4420424 f, uint3 x4unused4, Tuple_uint3_uint3_uint As)
 {
     Base3Broadcast344344544204247Tuple5Tuple5UInt327UInt327UInt3267Tuple5UInt327UInt327UInt3267UInt3266 x43;
     x43 = (Base3Broadcast344344544204247Tuple5Tuple5UInt327UInt327UInt3267Tuple5UInt327UInt327UInt3267UInt3266){f, As};
     Base3Broadcast344344544204247Tuple5Tuple5UInt327UInt327UInt3267Tuple5UInt327UInt327UInt3267UInt3266 _ssavalue_0;
     _ssavalue_0 = x43;
     Type5Val5366 _ssavalue_1;
     _ssavalue_1 = TYP_INST_Type5Val5366;
     return ntuple_184(_ssavalue_0, _ssavalue_1);
 }
 // (Base.Broadcast.broadcast_c, Tuple{##20#24,Type{Tuple},Tuple{UInt32,UInt32,UInt32},Tuple{UInt32,UInt32,UInt32},UInt32})
 uint3 broadcast_c_186(x4420424 f, Any x4unused4, uint3 A, Tuple_uint3_uint Bs)
 {
     Tuple_uint3_uint x44_apply_tmp41522;
     x44_apply_tmp41522 = Bs;
     Tuple_uint3_uint x44_apply_tmp41523;
     x44_apply_tmp41523 = Bs;
     return tuplebroadcast_185(f, first_tuple_177(A, (Tuple_uint3_uint){x44_apply_tmp41522.field1, x44_apply_tmp41522.field2}), (Tuple_uint3_uint3_uint){A, x44_apply_tmp41523.field1, x44_apply_tmp41523.field2});
 }
 // (broadcast, Tuple{##20#24,Tuple{UInt32,UInt32,UInt32},Tuple{UInt32,UInt32,UInt32},UInt32})
 uint3 broadcast_187(x4420424 f, uint3 A, Tuple_uint3_uint Bs)
 {
     Tuple_uint3_uint x44_apply_tmp41519;
     x44_apply_tmp41519 = Bs;
     Tuple_uint3_uint x44_apply_tmp41520;
     x44_apply_tmp41520 = Bs;
     return broadcast_c_186(f, containertype_177(A, x44_apply_tmp41519.field1, (uint){x44_apply_tmp41519.field2}), A, (Tuple_uint3_uint){x44_apply_tmp41520.field1, x44_apply_tmp41520.field2});
 }
 // Base.#getindex
 __constant int FUNC_INST_Base34getindex = 0;
 typedef int Base34getindex; // empty type emitted as an int
 // GPUArrays.#gpu_sub2ind
 __constant int FUNC_INST_GPUArrays34gpu_sub2ind = 0;
 typedef int GPUArrays34gpu_sub2ind; // empty type emitted as an int
 // GPUArrays.#_sub2ind
 __constant int FUNC_INST_GPUArrays34_sub2ind = 0;
 typedef int GPUArrays34_sub2ind; // empty type emitted as an int
 // (Base.argtail, Tuple{UInt32,UInt32})
 uint argtail_5(uint x, uint rest)
 {
     return rest;
 }
 // (Base.tail, Tuple{Tuple{UInt32,UInt32}})
 uint tail_6(uint2 x)
 {
     uint2 x44_apply_tmp41532;
     x44_apply_tmp41532 = x;
     return argtail_5(x44_apply_tmp41532.s0, (uint){x44_apply_tmp41532.s1});
 }
 // (GPUArrays._sub2ind, Tuple{Tuple{},UInt32,UInt32})
 uint _sub2ind_13(EmptyTuple_ x, uint L, uint ind)
 {
     return ind;
 }
 // (GPUArrays._sub2ind, Tuple{Tuple{UInt32},UInt32,UInt32,UInt32})
 uint _sub2ind_14(uint inds, uint L, uint ind, uint i, EmptyTuple_ I)
 {
     uint r1;
     r1 = inds;
     EmptyTuple_ x44_apply_tmp41533;
     x44_apply_tmp41533 = I;
     return _sub2ind_13(tail_12(inds), L * r1, ind + (i - (uint)(1)) * L);
 }
 // (GPUArrays._sub2ind, Tuple{Tuple{UInt32,UInt32},UInt32,UInt32,UInt32,UInt32})
 uint _sub2ind_15(uint2 inds, uint L, uint ind, uint i, uint I)
 {
     uint r1;
     r1 = inds.s0;
     uint x44_apply_tmp41531;
     x44_apply_tmp41531 = I;
     return _sub2ind_14(tail_6(inds), L * r1, ind + (i - (uint)(1)) * L, x44_apply_tmp41531, (EmptyTuple_){0.0f});
 }
 // (GPUArrays._sub2ind, Tuple{Tuple{UInt32,UInt32,UInt32},UInt32,UInt32,UInt32,UInt32,UInt32})
 uint _sub2ind_16(uint3 inds, uint L, uint ind, uint i, uint2 I)
 {
     uint r1;
     r1 = inds.s0;
     uint2 x44_apply_tmp41530;
     x44_apply_tmp41530 = I;
     return _sub2ind_15(tail_9(inds), L * r1, ind + (i - (uint)(1)) * L, x44_apply_tmp41530.s0, (uint){x44_apply_tmp41530.s1});
 }
 // (GPUArrays.gpu_sub2ind, Tuple{Tuple{UInt32,UInt32,UInt32},Tuple{UInt32,UInt32,UInt32}})
 uint gpu_sub2ind_17(uint3 dims, uint3 I)
 {
     uint3 x44_apply_tmp41529;
     x44_apply_tmp41529 = I;
     return _sub2ind_16((uint3)(dims), (uint)(1), (uint)(1), x44_apply_tmp41529.s0, (uint2){x44_apply_tmp41529.s1, x44_apply_tmp41529.s2});
 }
 // Base.#size
 __constant int FUNC_INST_Base34size = 0;
 typedef int Base34size; // empty type emitted as an int
 // (size, Tuple{CLArrays.DeviceArray{Tuple{Float32,Float32,Float32},3,Transpiler.CLIntrinsics.GlobalPointer{Tuple{Float32,Float32,Float32}}}})
 uint3 size_134(DeviceArray_float3_3___global1float3121 x)
 {
     return x.size;
 }
 // (map, Tuple{Type{UInt32},Tuple{UInt32,UInt32,UInt32}})
 uint3 map_89(Type5UInt326 f, uint3 t)
 {
     return (uint3){(uint)(t.s0), (uint)(t.s1), (uint)(t.s2)};
 }
 // (broadcast, Tuple{Type{UInt32},Tuple{UInt32,UInt32,UInt32}})
 uint3 broadcast_89(Type5UInt326 f, uint3 t, EmptyTuple_ ts)
 {
     EmptyTuple_ x44_apply_tmp41534;
     x44_apply_tmp41534 = ts;
     return map_89(f, t);
 }
 // Transpiler.#vload
 __constant int FUNC_INST_Transpiler34vload = 0;
 typedef int Transpiler34vload; // empty type emitted as an int
 // (Transpiler.vload, Tuple{Type{Tuple{Float32,Float32,Float32}},Transpiler.CLIntrinsics.GlobalPointer{Tuple{Float32,Float32,Float32}},UInt32})
 float3 vload_35(Any x4unused4, __global float3 *  a, uint i)
 {
     return (float3)(vload3(i - (uint)(1), (__global float * )(a)));
 }
 // (getindex, Tuple{Transpiler.CLIntrinsics.GlobalPointer{Tuple{Float32,Float32,Float32}},UInt32})
 float3 getindex_36(__global float3 *  a, uint i)
 {
     return vload_35(TYP_INST_Type5Tuple5Float327Float327Float3266, a, i);
 }
 // (getindex, Tuple{CLArrays.DeviceArray{Tuple{Float32,Float32,Float32},3,Transpiler.CLIntrinsics.GlobalPointer{Tuple{Float32,Float32,Float32}}},UInt32,UInt32,UInt32})
 float3 getindex_135(DeviceArray_float3_3___global1float3121 x, uint3 i)
 {
     uint ilin;
     ilin = gpu_sub2ind_17(size_134(x), broadcast_89(TYP_INST_Type5UInt326, i, (EmptyTuple_){0.0f}));
     return getindex_36(x.ptr, ilin);
 }
 // Type{##21#25}
 typedef int Type544214256; // placeholder type instance
 __constant Type544214256 TYP_INST_Type544214256 = 0;

 // Tuple{Tuple{UInt32,UInt32,UInt32},Tuple{Float32,Float32,Float32}}
 struct  __attribute__ ((packed)) TYPTuple_uint3_float3{
     uint3 field1;
     float3 field2;
 };
 typedef struct TYPTuple_uint3_float3 Tuple_uint3_float3;

 // Tuple{Tuple{Float32,Float32,Float32}}
 struct  __attribute__ ((packed)) TYPTuple_float3{
     float3 field1;
 };
 typedef struct TYPTuple_float3 Tuple_float3;

 // Tuple{Float32,UInt32,Float32}
 struct  __attribute__ ((packed)) TYPTuple_float_uint_float{
     float field1;
     uint field2;
     float field3;
 };
 typedef struct TYPTuple_float_uint_float Tuple_float_uint_float;

 // Base.#heads
 __constant int FUNC_INST_Base34heads = 0;
 typedef int Base34heads; // empty type emitted as an int
 // Tuple{Tuple{Float32,Float32,Float32},Tuple{UInt32,UInt32,UInt32},Tuple{Float32,Float32,Float32}}
 struct  __attribute__ ((packed)) TYPTuple_float3_uint3_float3{
     float3 field1;
     uint3 field2;
     float3 field3;
 };
 typedef struct TYPTuple_float3_uint3_float3 Tuple_float3_uint3_float3;

 // Base.##28#29
 __constant int FUNC_INST_Base34428429 = 0;
 typedef int Base34428429; // empty type emitted as an int
 // Type{Base.##28#29}
 typedef int Type5Base344284296; // placeholder type instance
 __constant Type5Base344284296 TYP_INST_Type5Base344284296 = 0;

 // (Base.#28, Tuple{Tuple{Float32,Float32,Float32}})
 float x428_122(float3 t)
 {
     return t.s0;
 }
 // (Base.#28, Tuple{Tuple{UInt32,UInt32,UInt32}})
 uint x428_9(uint3 t)
 {
     return t.s0;
 }
 // (map, Tuple{Base.##28#29,Tuple{Tuple{Float32,Float32,Float32},Tuple{UInt32,UInt32,UInt32},Tuple{Float32,Float32,Float32}}})
 Tuple_float_uint_float map_188(Base34428429 f, Tuple_float3_uint3_float3 t)
 {
     return (Tuple_float_uint_float){x428_122(t.field1), x428_9(t.field2), x428_122(t.field3)};
 }
 // (Base.heads, Tuple{Tuple{Float32,Float32,Float32},Tuple{UInt32,UInt32,UInt32},Tuple{Float32,Float32,Float32}})
 Tuple_float_uint_float heads_189(Tuple_float3_uint3_float3 ts)
 {
     Base34428429 x428;
     x428 = (Base34428429){0.0f};
     Base34428429 _ssavalue_0;
     _ssavalue_0 = x428;
     return map_188(_ssavalue_0, ts);
 }
 // (#21, Tuple{Float32,UInt32,Float32})
 float x421_190(float x4temp4, uint x4temp41, float x4temp42)
 {
     float x4temp44;
     float x4temp43;
     x4temp43 = x4temp41 - 1.0f;
     x4temp44 = x4temp43 * x4temp42;
     return x4temp4 - x4temp44;
 }
 // Tuple{Tuple{Float32,Float32},Tuple{UInt32,UInt32},Tuple{Float32,Float32}}
 struct  __attribute__ ((packed)) TYPTuple_float2_uint2_float2{
     float2 field1;
     uint2 field2;
     float2 field3;
 };
 typedef struct TYPTuple_float2_uint2_float2 Tuple_float2_uint2_float2;

 // Base.#tails
 __constant int FUNC_INST_Base34tails = 0;
 typedef int Base34tails; // empty type emitted as an int
 // (map, Tuple{Base.#tail,Tuple{Tuple{Float32,Float32,Float32},Tuple{UInt32,UInt32,UInt32},Tuple{Float32,Float32,Float32}}})
 Tuple_float2_uint2_float2 map_191(Base34tail f, Tuple_float3_uint3_float3 t)
 {
     return (Tuple_float2_uint2_float2){tail_122(t.field1), tail_9(t.field2), tail_122(t.field3)};
 }
 // (Base.tails, Tuple{Tuple{Float32,Float32,Float32},Tuple{UInt32,UInt32,UInt32},Tuple{Float32,Float32,Float32}})
 Tuple_float2_uint2_float2 tails_189(Tuple_float3_uint3_float3 ts)
 {
     return map_191(FUNC_INST_Base34tail, ts);
 }
 // Tuple{Tuple{Float32,Float32}}
 struct  __attribute__ ((packed)) TYPTuple_float2{
     float2 field1;
 };
 typedef struct TYPTuple_float2 Tuple_float2;

 // (Base.#28, Tuple{Tuple{Float32,Float32}})
 float x428_57(float2 t)
 {
     return t.s0;
 }
 // (Base.#28, Tuple{Tuple{UInt32,UInt32}})
 uint x428_6(uint2 t)
 {
     return t.s0;
 }
 // (map, Tuple{Base.##28#29,Tuple{Tuple{Float32,Float32},Tuple{UInt32,UInt32},Tuple{Float32,Float32}}})
 Tuple_float_uint_float map_192(Base34428429 f, Tuple_float2_uint2_float2 t)
 {
     return (Tuple_float_uint_float){x428_57(t.field1), x428_6(t.field2), x428_57(t.field3)};
 }
 // (Base.heads, Tuple{Tuple{Float32,Float32},Tuple{UInt32,UInt32},Tuple{Float32,Float32}})
 Tuple_float_uint_float heads_193(Tuple_float2_uint2_float2 ts)
 {
     Base34428429 x428;
     x428 = (Base34428429){0.0f};
     Base34428429 _ssavalue_0;
     _ssavalue_0 = x428;
     return map_192(_ssavalue_0, ts);
 }
 // Tuple{Tuple{Float32},Tuple{UInt32},Tuple{Float32}}
 struct  __attribute__ ((packed)) TYPTuple_float_uint_float{
     float field1;
     uint field2;
     float field3;
 };
 typedef struct TYPTuple_float_uint_float Tuple_float_uint_float;

 // (Base.argtail, Tuple{Float32,Float32})
 float argtail_24(float x, float rest)
 {
     return rest;
 }
 // (Base.tail, Tuple{Tuple{Float32,Float32}})
 float tail_57(float2 x)
 {
     float2 x44_apply_tmp41546;
     x44_apply_tmp41546 = x;
     return argtail_24(x44_apply_tmp41546.s0, (float){x44_apply_tmp41546.s1});
 }
 // (map, Tuple{Base.#tail,Tuple{Tuple{Float32,Float32},Tuple{UInt32,UInt32},Tuple{Float32,Float32}}})
 Tuple_float_uint_float map_194(Base34tail f, Tuple_float2_uint2_float2 t)
 {
     return (Tuple_float_uint_float){tail_57(t.field1), tail_6(t.field2), tail_57(t.field3)};
 }
 // (Base.tails, Tuple{Tuple{Float32,Float32},Tuple{UInt32,UInt32},Tuple{Float32,Float32}})
 Tuple_float_uint_float tails_193(Tuple_float2_uint2_float2 ts)
 {
     return map_194(FUNC_INST_Base34tail, ts);
 }
 // Tuple{Tuple{Float32}}
 struct  __attribute__ ((packed)) TYPTuple_float{
     float field1;
 };
 typedef struct TYPTuple_float Tuple_float;

 // (Base.#28, Tuple{Tuple{Float32}})
 float x428_107(float t)
 {
     return t;
 }
 // (Base.#28, Tuple{Tuple{UInt32}})
 uint x428_12(uint t)
 {
     return t;
 }
 // (map, Tuple{Base.##28#29,Tuple{Tuple{Float32},Tuple{UInt32},Tuple{Float32}}})
 Tuple_float_uint_float map_195(Base34428429 f, Tuple_float_uint_float t)
 {
     return (Tuple_float_uint_float){x428_107(t.field1), x428_12(t.field2), x428_107(t.field3)};
 }
 // (Base.heads, Tuple{Tuple{Float32},Tuple{UInt32},Tuple{Float32}})
 Tuple_float_uint_float heads_196(Tuple_float_uint_float ts)
 {
     Base34428429 x428;
     x428 = (Base34428429){0.0f};
     Base34428429 _ssavalue_0;
     _ssavalue_0 = x428;
     return map_195(_ssavalue_0, ts);
 }
 // Tuple{Tuple{},Tuple{},Tuple{}}
 struct  __attribute__ ((packed)) TYPTuple_EmptyTuple__EmptyTuple__EmptyTuple_{
     EmptyTuple_ field1;
     EmptyTuple_ field2;
     EmptyTuple_ field3;
 };
 typedef struct TYPTuple_EmptyTuple__EmptyTuple__EmptyTuple_ Tuple_EmptyTuple__EmptyTuple__EmptyTuple_;

 // (Base.argtail, Tuple{Float32})
 EmptyTuple_ argtail_22(float x, EmptyTuple_ rest)
 {
     return rest;
 }
 // (Base.tail, Tuple{Tuple{Float32}})
 EmptyTuple_ tail_107(float x)
 {
     float x44_apply_tmp41552;
     x44_apply_tmp41552 = x;
     return argtail_22(x44_apply_tmp41552, (EmptyTuple_){0.0f});
 }
 // (map, Tuple{Base.#tail,Tuple{Tuple{Float32},Tuple{UInt32},Tuple{Float32}}})
 Tuple_EmptyTuple__EmptyTuple__EmptyTuple_ map_197(Base34tail f, Tuple_float_uint_float t)
 {
     return (Tuple_EmptyTuple__EmptyTuple__EmptyTuple_){tail_107(t.field1), tail_12(t.field2), tail_107(t.field3)};
 }
 // (Base.tails, Tuple{Tuple{Float32},Tuple{UInt32},Tuple{Float32}})
 Tuple_EmptyTuple__EmptyTuple__EmptyTuple_ tails_196(Tuple_float_uint_float ts)
 {
     return map_197(FUNC_INST_Base34tail, ts);
 }
 // (map, Tuple{##21#25,Tuple{},Tuple{},Tuple{}})
 EmptyTuple_ map_198(x4421425 f, Tuple_EmptyTuple__EmptyTuple__EmptyTuple_ x4unused4)
 {
     return (EmptyTuple_){0.0f};
 }
 // (map, Tuple{##21#25,Tuple{Float32},Tuple{UInt32},Tuple{Float32}})
 float map_199(x4421425 f, float t1, uint t2, Tuple_float ts)
 {
     Tuple_float_uint_float x44_apply_tmp41548;
     x44_apply_tmp41548 = Sugar.InlineNode(Any[:(_6::Tuple{Tuple{Float32}}::Tuple{Tuple{Float32}}), :(_6::Tuple{Tuple{Float32}} = _5::Tuple{Tuple{Float32}})], :((Base.heads)((Tuple{Tuple{Float32},Tuple{UInt32},Tuple{Float32}}){_3::Tuple{Float32}, _4::Tuple{UInt32}, (getfield)(_6::Tuple{Tuple{Float32}}, field1)::Tuple{Float32}}::Tuple{Tuple{Float32},Tuple{UInt32},Tuple{Float32}})::Tuple{Float32,UInt32,Float32}));
     EmptyTuple_ x44_apply_tmp41551;
     x44_apply_tmp41551 = Sugar.InlineNode(Any[:(_9::Tuple{Tuple{},Tuple{},Tuple{}}::Tuple{Tuple{},Tuple{},Tuple{}}), :(_9::Tuple{Tuple{},Tuple{},Tuple{}} = Sugar.InlineNode(Any[:(_8::Tuple{Tuple{Float32}}::Tuple{Tuple{Float32}}), :(_8::Tuple{Tuple{Float32}} = _5::Tuple{Tuple{Float32}})], :((Base.tails)((Tuple{Tuple{Float32},Tuple{UInt32},Tuple{Float32}}){_3::Tuple{Float32}, _4::Tuple{UInt32}, (getfield)(_8::Tuple{Tuple{Float32}}, field1)::Tuple{Float32}}::Tuple{Tuple{Float32},Tuple{UInt32},Tuple{Float32}})::Tuple{Tuple{},Tuple{},Tuple{}})))], :((map)(_2::##21#25, (Tuple{Tuple{},Tuple{},Tuple{}}){(getfield)(_9::Tuple{Tuple{},Tuple{},Tuple{}}, field1)::Tuple{}, (getfield)(_9::Tuple{Tuple{},Tuple{},Tuple{}}, field2)::Tuple{}, (getfield)(_9::Tuple{Tuple{},Tuple{},Tuple{}}, field3)::Tuple{}}::Tuple{Tuple{},Tuple{},Tuple{}})::Tuple{}));
     return (float){x421_190(x44_apply_tmp41548.field1, x44_apply_tmp41548.field2, x44_apply_tmp41548.field3)};
 }
 // (map, Tuple{##21#25,Tuple{Float32,Float32},Tuple{UInt32,UInt32},Tuple{Float32,Float32}})
 float2 map_200(x4421425 f, float2 t1, uint2 t2, Tuple_float2 ts)
 {
     Tuple_float_uint_float x44_apply_tmp41542;
     x44_apply_tmp41542 = Sugar.InlineNode(Any[:(_6::Tuple{Tuple{Float32,Float32}}::Tuple{Tuple{Float32,Float32}}), :(_6::Tuple{Tuple{Float32,Float32}} = _5::Tuple{Tuple{Float32,Float32}})], :((Base.heads)((Tuple{Tuple{Float32,Float32},Tuple{UInt32,UInt32},Tuple{Float32,Float32}}){_3::Tuple{Float32,Float32}, _4::Tuple{UInt32,UInt32}, (getfield)(_6::Tuple{Tuple{Float32,Float32}}, field1)::Tuple{Float32,Float32}}::Tuple{Tuple{Float32,Float32},Tuple{UInt32,UInt32},Tuple{Float32,Float32}})::Tuple{Float32,UInt32,Float32}));
     float x44_apply_tmp41545;
     x44_apply_tmp41545 = Sugar.InlineNode(Any[:(_9::Tuple{Tuple{Float32},Tuple{UInt32},Tuple{Float32}}::Tuple{Tuple{Float32},Tuple{UInt32},Tuple{Float32}}), :(_9::Tuple{Tuple{Float32},Tuple{UInt32},Tuple{Float32}} = Sugar.InlineNode(Any[:(_8::Tuple{Tuple{Float32,Float32}}::Tuple{Tuple{Float32,Float32}}), :(_8::Tuple{Tuple{Float32,Float32}} = _5::Tuple{Tuple{Float32,Float32}})], :((Base.tails)((Tuple{Tuple{Float32,Float32},Tuple{UInt32,UInt32},Tuple{Float32,Float32}}){_3::Tuple{Float32,Float32}, _4::Tuple{UInt32,UInt32}, (getfield)(_8::Tuple{Tuple{Float32,Float32}}, field1)::Tuple{Float32,Float32}}::Tuple{Tuple{Float32,Float32},Tuple{UInt32,UInt32},Tuple{Float32,Float32}})::Tuple{Tuple{Float32},Tuple{UInt32},Tuple{Float32}})))], :((map)(_2::##21#25, (getfield)(_9::Tuple{Tuple{Float32},Tuple{UInt32},Tuple{Float32}}, field1)::Tuple{Float32}, (getfield)(_9::Tuple{Tuple{Float32},Tuple{UInt32},Tuple{Float32}}, field2)::Tuple{UInt32}, (Tuple{Tuple{Float32}}){(getfield)(_9::Tuple{Tuple{Float32},Tuple{UInt32},Tuple{Float32}}, field3)::Tuple{Float32}}::Tuple{Tuple{Float32}})::Tuple{Float32}));
     return (float2){x421_190(x44_apply_tmp41542.field1, x44_apply_tmp41542.field2, x44_apply_tmp41542.field3), x44_apply_tmp41545};
 }
 // (map, Tuple{##21#25,Tuple{Float32,Float32,Float32},Tuple{UInt32,UInt32,UInt32},Tuple{Float32,Float32,Float32}})
 float3 map_201(x4421425 f, float3 t1, uint3 t2, Tuple_float3 ts)
 {
     Tuple_float_uint_float x44_apply_tmp41537;
     x44_apply_tmp41537 = Sugar.InlineNode(Any[:(_6::Tuple{Tuple{Float32,Float32,Float32}}::Tuple{Tuple{Float32,Float32,Float32}}), :(_6::Tuple{Tuple{Float32,Float32,Float32}} = _5::Tuple{Tuple{Float32,Float32,Float32}})], :((Base.heads)((Tuple{Tuple{Float32,Float32,Float32},Tuple{UInt32,UInt32,UInt32},Tuple{Float32,Float32,Float32}}){_3::Tuple{Float32,Float32,Float32}, _4::Tuple{UInt32,UInt32,UInt32}, (getfield)(_6::Tuple{Tuple{Float32,Float32,Float32}}, field1)::Tuple{Float32,Float32,Float32}}::Tuple{Tuple{Float32,Float32,Float32},Tuple{UInt32,UInt32,UInt32},Tuple{Float32,Float32,Float32}})::Tuple{Float32,UInt32,Float32}));
     float2 x44_apply_tmp41540;
     x44_apply_tmp41540 = Sugar.InlineNode(Any[:(_9::Tuple{Tuple{Float32,Float32},Tuple{UInt32,UInt32},Tuple{Float32,Float32}}::Tuple{Tuple{Float32,Float32},Tuple{UInt32,UInt32},Tuple{Float32,Float32}}), :(_9::Tuple{Tuple{Float32,Float32},Tuple{UInt32,UInt32},Tuple{Float32,Float32}} = Sugar.InlineNode(Any[:(_8::Tuple{Tuple{Float32,Float32,Float32}}::Tuple{Tuple{Float32,Float32,Float32}}), :(_8::Tuple{Tuple{Float32,Float32,Float32}} = _5::Tuple{Tuple{Float32,Float32,Float32}})], :((Base.tails)((Tuple{Tuple{Float32,Float32,Float32},Tuple{UInt32,UInt32,UInt32},Tuple{Float32,Float32,Float32}}){_3::Tuple{Float32,Float32,Float32}, _4::Tuple{UInt32,UInt32,UInt32}, (getfield)(_8::Tuple{Tuple{Float32,Float32,Float32}}, field1)::Tuple{Float32,Float32,Float32}}::Tuple{Tuple{Float32,Float32,Float32},Tuple{UInt32,UInt32,UInt32},Tuple{Float32,Float32,Float32}})::Tuple{Tuple{Float32,Float32},Tuple{UInt32,UInt32},Tuple{Float32,Float32}})))], :((map)(_2::##21#25, (getfield)(_9::Tuple{Tuple{Float32,Float32},Tuple{UInt32,UInt32},Tuple{Float32,Float32}}, field1)::Tuple{Float32,Float32}, (getfield)(_9::Tuple{Tuple{Float32,Float32},Tuple{UInt32,UInt32},Tuple{Float32,Float32}}, field2)::Tuple{UInt32,UInt32}, (Tuple{Tuple{Float32,Float32}}){(getfield)(_9::Tuple{Tuple{Float32,Float32},Tuple{UInt32,UInt32},Tuple{Float32,Float32}}, field3)::Tuple{Float32,Float32}}::Tuple{Tuple{Float32,Float32}})::Tuple{Float32,Float32}));
     return (float3){x421_190(x44_apply_tmp41537.field1, x44_apply_tmp41537.field2, x44_apply_tmp41537.field3), x44_apply_tmp41540.s0, x44_apply_tmp41540.s1};
 }
 // (broadcast, Tuple{##21#25,Tuple{Float32,Float32,Float32},Tuple{UInt32,UInt32,UInt32},Tuple{Float32,Float32,Float32}})
 float3 broadcast_201(x4421425 f, float3 t, Tuple_uint3_float3 ts)
 {
     Tuple_uint3_float3 x44_apply_tmp41535;
     x44_apply_tmp41535 = ts;
     return map_201(f, t, x44_apply_tmp41535.field1, (Tuple_float3){x44_apply_tmp41535.field2});
 }
 // Type{##22#26}
 typedef int Type544224266; // placeholder type instance
 __constant Type544224266 TYP_INST_Type544224266 = 0;

 // NTuple{5,Tuple{Float32,Float32,Float32}}
 struct  __attribute__ ((packed)) TYPTuple_float3_float3_float3_float3_float3{
     float3 field1;
     float3 field2;
     float3 field3;
     float3 field4;
     float3 field5;
 };
 typedef struct TYPTuple_float3_float3_float3_float3_float3 Tuple_float3_float3_float3_float3_float3;

 // NTuple{4,Tuple{Float32,Float32,Float32}}
 struct  __attribute__ ((packed)) TYPTuple_float3_float3_float3_float3{
     float3 field1;
     float3 field2;
     float3 field3;
     float3 field4;
 };
 typedef struct TYPTuple_float3_float3_float3_float3 Tuple_float3_float3_float3_float3;

 // NTuple{6,Float32}
 struct  __attribute__ ((packed)) TYPTuple_float_float_float_float_float_float{
     float field1;
     float field2;
     float field3;
     float field4;
     float field5;
     float field6;
 };
 typedef struct TYPTuple_float_float_float_float_float_float Tuple_float_float_float_float_float_float;

 // NTuple{6,Tuple{Float32,Float32,Float32}}
 struct  __attribute__ ((packed)) TYPTuple_float3_float3_float3_float3_float3_float3{
     float3 field1;
     float3 field2;
     float3 field3;
     float3 field4;
     float3 field5;
     float3 field6;
 };
 typedef struct TYPTuple_float3_float3_float3_float3_float3_float3 Tuple_float3_float3_float3_float3_float3_float3;

 // NTuple{5,Float32}
 struct  __attribute__ ((packed)) TYPTuple_float_float_float_float_float{
     float field1;
     float field2;
     float field3;
     float field4;
     float field5;
 };
 typedef struct TYPTuple_float_float_float_float_float Tuple_float_float_float_float_float;

 // Tuple{Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32}}
 struct  __attribute__ ((packed)) TYPTuple_float3_float3_float3{
     float3 field1;
     float3 field2;
     float3 field3;
 };
 typedef struct TYPTuple_float3_float3_float3 Tuple_float3_float3_float3;

 // (map, Tuple{Base.##28#29,Tuple{Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32}}})
 float3 map_202(Base34428429 f, Tuple_float3_float3_float3 t)
 {
     return (float3){x428_122(t.field1), x428_122(t.field2), x428_122(t.field3)};
 }
 // (Base.argtail, NTuple{4,Tuple{Float32,Float32,Float32}})
 Tuple_float3_float3_float3 argtail_203(float3 x, Tuple_float3_float3_float3 rest)
 {
     return rest;
 }
 // (Base.tail, Tuple{NTuple{4,Tuple{Float32,Float32,Float32}}})
 Tuple_float3_float3_float3 tail_204(Tuple_float3_float3_float3_float3 x)
 {
     Tuple_float3_float3_float3_float3 x44_apply_tmp41562;
     x44_apply_tmp41562 = x;
     return argtail_203(x44_apply_tmp41562.field1, (Tuple_float3_float3_float3){x44_apply_tmp41562.field2, x44_apply_tmp41562.field3, x44_apply_tmp41562.field4});
 }
 // (map, Tuple{Base.##28#29,NTuple{4,Tuple{Float32,Float32,Float32}}})
 float4 map_205(Base34428429 f, Tuple_float3_float3_float3_float3 t)
 {
     float3 x44_apply_tmp41561;
     x44_apply_tmp41561 = map_202(f, tail_204(t));
     return (float4){x428_122(t.field1), x44_apply_tmp41561.s0, x44_apply_tmp41561.s1, x44_apply_tmp41561.s2};
 }
 // (Base.argtail, NTuple{5,Tuple{Float32,Float32,Float32}})
 Tuple_float3_float3_float3_float3 argtail_206(float3 x, Tuple_float3_float3_float3_float3 rest)
 {
     return rest;
 }
 // (Base.tail, Tuple{NTuple{5,Tuple{Float32,Float32,Float32}}})
 Tuple_float3_float3_float3_float3 tail_207(Tuple_float3_float3_float3_float3_float3 x)
 {
     Tuple_float3_float3_float3_float3_float3 x44_apply_tmp41563;
     x44_apply_tmp41563 = x;
     return argtail_206(x44_apply_tmp41563.field1, (Tuple_float3_float3_float3_float3){x44_apply_tmp41563.field2, x44_apply_tmp41563.field3, x44_apply_tmp41563.field4, x44_apply_tmp41563.field5});
 }
 // (map, Tuple{Base.##28#29,NTuple{5,Tuple{Float32,Float32,Float32}}})
 Tuple_float_float_float_float_float map_208(Base34428429 f, Tuple_float3_float3_float3_float3_float3 t)
 {
     float4 x44_apply_tmp41560;
     x44_apply_tmp41560 = map_205(f, tail_207(t));
     return (Tuple_float_float_float_float_float){x428_122(t.field1), x44_apply_tmp41560.s0, x44_apply_tmp41560.s1, x44_apply_tmp41560.s2, x44_apply_tmp41560.s3};
 }
 // (Base.argtail, NTuple{6,Tuple{Float32,Float32,Float32}})
 Tuple_float3_float3_float3_float3_float3 argtail_209(float3 x, Tuple_float3_float3_float3_float3_float3 rest)
 {
     return rest;
 }
 // (Base.tail, Tuple{NTuple{6,Tuple{Float32,Float32,Float32}}})
 Tuple_float3_float3_float3_float3_float3 tail_210(Tuple_float3_float3_float3_float3_float3_float3 x)
 {
     Tuple_float3_float3_float3_float3_float3_float3 x44_apply_tmp41564;
     x44_apply_tmp41564 = x;
     return argtail_209(x44_apply_tmp41564.field1, (Tuple_float3_float3_float3_float3_float3){x44_apply_tmp41564.field2, x44_apply_tmp41564.field3, x44_apply_tmp41564.field4, x44_apply_tmp41564.field5, x44_apply_tmp41564.field6});
 }
 // (map, Tuple{Base.##28#29,NTuple{6,Tuple{Float32,Float32,Float32}}})
 Tuple_float_float_float_float_float_float map_211(Base34428429 f, Tuple_float3_float3_float3_float3_float3_float3 t)
 {
     Tuple_float_float_float_float_float x44_apply_tmp41559;
     x44_apply_tmp41559 = map_208(f, tail_210(t));
     return (Tuple_float_float_float_float_float_float){x428_122(t.field1), x44_apply_tmp41559.field1, x44_apply_tmp41559.field2, x44_apply_tmp41559.field3, x44_apply_tmp41559.field4, x44_apply_tmp41559.field5};
 }
 // (Base.heads, NTuple{6,Tuple{Float32,Float32,Float32}})
 Tuple_float_float_float_float_float_float heads_209(Tuple_float3_float3_float3_float3_float3_float3 ts)
 {
     Base34428429 x428;
     x428 = (Base34428429){0.0f};
     Base34428429 _ssavalue_0;
     _ssavalue_0 = x428;
     return map_211(_ssavalue_0, ts);
 }
 // (#22, NTuple{6,Float32})
 float x422_212(float x4temp4, float x4temp41, float x4temp42, float x4temp43, float x4temp44, float x4temp45)
 {
     float x4temp416;
     float x4temp415;
     float x4temp414;
     float x4temp413;
     float x4temp412;
     float x4temp411;
     float x4temp410;
     float x4temp49;
     float x4temp48;
     float x4temp47;
     float x4temp46;
     x4temp410 = 1.0f - x4temp4;
     x4temp46 = 1.0f - x4temp41;
     x4temp48 = x4temp46 * x4temp42;
     x4temp47 = x4temp41 * x4temp43;
     x4temp49 = x4temp48 + x4temp47;
     x4temp416 = x4temp410 * x4temp49;
     x4temp411 = 1.0f - x4temp41;
     x4temp413 = x4temp411 * x4temp44;
     x4temp412 = x4temp41 * x4temp45;
     x4temp414 = x4temp413 + x4temp412;
     x4temp415 = x4temp4 * x4temp414;
     return x4temp416 + x4temp415;
 }
 // NTuple{6,Tuple{Float32,Float32}}
 struct  __attribute__ ((packed)) TYPTuple_float2_float2_float2_float2_float2_float2{
     float2 field1;
     float2 field2;
     float2 field3;
     float2 field4;
     float2 field5;
     float2 field6;
 };
 typedef struct TYPTuple_float2_float2_float2_float2_float2_float2 Tuple_float2_float2_float2_float2_float2_float2;

 // NTuple{5,Tuple{Float32,Float32}}
 struct  __attribute__ ((packed)) TYPTuple_float2_float2_float2_float2_float2{
     float2 field1;
     float2 field2;
     float2 field3;
     float2 field4;
     float2 field5;
 };
 typedef struct TYPTuple_float2_float2_float2_float2_float2 Tuple_float2_float2_float2_float2_float2;

 // NTuple{4,Tuple{Float32,Float32}}
 struct  __attribute__ ((packed)) TYPTuple_float2_float2_float2_float2{
     float2 field1;
     float2 field2;
     float2 field3;
     float2 field4;
 };
 typedef struct TYPTuple_float2_float2_float2_float2 Tuple_float2_float2_float2_float2;

 // Tuple{Tuple{Float32,Float32},Tuple{Float32,Float32},Tuple{Float32,Float32}}
 struct  __attribute__ ((packed)) TYPTuple_float2_float2_float2{
     float2 field1;
     float2 field2;
     float2 field3;
 };
 typedef struct TYPTuple_float2_float2_float2 Tuple_float2_float2_float2;

 // (map, Tuple{Base.#tail,Tuple{Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32}}})
 Tuple_float2_float2_float2 map_213(Base34tail f, Tuple_float3_float3_float3 t)
 {
     return (Tuple_float2_float2_float2){tail_122(t.field1), tail_122(t.field2), tail_122(t.field3)};
 }
 // (map, Tuple{Base.#tail,NTuple{4,Tuple{Float32,Float32,Float32}}})
 Tuple_float2_float2_float2_float2 map_214(Base34tail f, Tuple_float3_float3_float3_float3 t)
 {
     Tuple_float2_float2_float2 x44_apply_tmp41567;
     x44_apply_tmp41567 = map_213(f, tail_204(t));
     return (Tuple_float2_float2_float2_float2){tail_122(t.field1), x44_apply_tmp41567.field1, x44_apply_tmp41567.field2, x44_apply_tmp41567.field3};
 }
 // (map, Tuple{Base.#tail,NTuple{5,Tuple{Float32,Float32,Float32}}})
 Tuple_float2_float2_float2_float2_float2 map_215(Base34tail f, Tuple_float3_float3_float3_float3_float3 t)
 {
     Tuple_float2_float2_float2_float2 x44_apply_tmp41566;
     x44_apply_tmp41566 = map_214(f, tail_207(t));
     return (Tuple_float2_float2_float2_float2_float2){tail_122(t.field1), x44_apply_tmp41566.field1, x44_apply_tmp41566.field2, x44_apply_tmp41566.field3, x44_apply_tmp41566.field4};
 }
 // (map, Tuple{Base.#tail,NTuple{6,Tuple{Float32,Float32,Float32}}})
 Tuple_float2_float2_float2_float2_float2_float2 map_216(Base34tail f, Tuple_float3_float3_float3_float3_float3_float3 t)
 {
     Tuple_float2_float2_float2_float2_float2 x44_apply_tmp41565;
     x44_apply_tmp41565 = map_215(f, tail_210(t));
     return (Tuple_float2_float2_float2_float2_float2_float2){tail_122(t.field1), x44_apply_tmp41565.field1, x44_apply_tmp41565.field2, x44_apply_tmp41565.field3, x44_apply_tmp41565.field4, x44_apply_tmp41565.field5};
 }
 // (Base.tails, NTuple{6,Tuple{Float32,Float32,Float32}})
 Tuple_float2_float2_float2_float2_float2_float2 tails_209(Tuple_float3_float3_float3_float3_float3_float3 ts)
 {
     return map_216(FUNC_INST_Base34tail, ts);
 }
 // (map, Tuple{Base.##28#29,Tuple{Tuple{Float32,Float32},Tuple{Float32,Float32},Tuple{Float32,Float32}}})
 float3 map_217(Base34428429 f, Tuple_float2_float2_float2 t)
 {
     return (float3){x428_57(t.field1), x428_57(t.field2), x428_57(t.field3)};
 }
 // (Base.argtail, NTuple{4,Tuple{Float32,Float32}})
 Tuple_float2_float2_float2 argtail_218(float2 x, Tuple_float2_float2_float2 rest)
 {
     return rest;
 }
 // (Base.tail, Tuple{NTuple{4,Tuple{Float32,Float32}}})
 Tuple_float2_float2_float2 tail_219(Tuple_float2_float2_float2_float2 x)
 {
     Tuple_float2_float2_float2_float2 x44_apply_tmp41576;
     x44_apply_tmp41576 = x;
     return argtail_218(x44_apply_tmp41576.field1, (Tuple_float2_float2_float2){x44_apply_tmp41576.field2, x44_apply_tmp41576.field3, x44_apply_tmp41576.field4});
 }
 // (map, Tuple{Base.##28#29,NTuple{4,Tuple{Float32,Float32}}})
 float4 map_220(Base34428429 f, Tuple_float2_float2_float2_float2 t)
 {
     float3 x44_apply_tmp41575;
     x44_apply_tmp41575 = map_217(f, tail_219(t));
     return (float4){x428_57(t.field1), x44_apply_tmp41575.s0, x44_apply_tmp41575.s1, x44_apply_tmp41575.s2};
 }
 // (Base.argtail, NTuple{5,Tuple{Float32,Float32}})
 Tuple_float2_float2_float2_float2 argtail_221(float2 x, Tuple_float2_float2_float2_float2 rest)
 {
     return rest;
 }
 // (Base.tail, Tuple{NTuple{5,Tuple{Float32,Float32}}})
 Tuple_float2_float2_float2_float2 tail_222(Tuple_float2_float2_float2_float2_float2 x)
 {
     Tuple_float2_float2_float2_float2_float2 x44_apply_tmp41577;
     x44_apply_tmp41577 = x;
     return argtail_221(x44_apply_tmp41577.field1, (Tuple_float2_float2_float2_float2){x44_apply_tmp41577.field2, x44_apply_tmp41577.field3, x44_apply_tmp41577.field4, x44_apply_tmp41577.field5});
 }
 // (map, Tuple{Base.##28#29,NTuple{5,Tuple{Float32,Float32}}})
 Tuple_float_float_float_float_float map_223(Base34428429 f, Tuple_float2_float2_float2_float2_float2 t)
 {
     float4 x44_apply_tmp41574;
     x44_apply_tmp41574 = map_220(f, tail_222(t));
     return (Tuple_float_float_float_float_float){x428_57(t.field1), x44_apply_tmp41574.s0, x44_apply_tmp41574.s1, x44_apply_tmp41574.s2, x44_apply_tmp41574.s3};
 }
 // (Base.argtail, NTuple{6,Tuple{Float32,Float32}})
 Tuple_float2_float2_float2_float2_float2 argtail_224(float2 x, Tuple_float2_float2_float2_float2_float2 rest)
 {
     return rest;
 }
 // (Base.tail, Tuple{NTuple{6,Tuple{Float32,Float32}}})
 Tuple_float2_float2_float2_float2_float2 tail_225(Tuple_float2_float2_float2_float2_float2_float2 x)
 {
     Tuple_float2_float2_float2_float2_float2_float2 x44_apply_tmp41578;
     x44_apply_tmp41578 = x;
     return argtail_224(x44_apply_tmp41578.field1, (Tuple_float2_float2_float2_float2_float2){x44_apply_tmp41578.field2, x44_apply_tmp41578.field3, x44_apply_tmp41578.field4, x44_apply_tmp41578.field5, x44_apply_tmp41578.field6});
 }
 // (map, Tuple{Base.##28#29,NTuple{6,Tuple{Float32,Float32}}})
 Tuple_float_float_float_float_float_float map_226(Base34428429 f, Tuple_float2_float2_float2_float2_float2_float2 t)
 {
     Tuple_float_float_float_float_float x44_apply_tmp41573;
     x44_apply_tmp41573 = map_223(f, tail_225(t));
     return (Tuple_float_float_float_float_float_float){x428_57(t.field1), x44_apply_tmp41573.field1, x44_apply_tmp41573.field2, x44_apply_tmp41573.field3, x44_apply_tmp41573.field4, x44_apply_tmp41573.field5};
 }
 // (Base.heads, NTuple{6,Tuple{Float32,Float32}})
 Tuple_float_float_float_float_float_float heads_224(Tuple_float2_float2_float2_float2_float2_float2 ts)
 {
     Base34428429 x428;
     x428 = (Base34428429){0.0f};
     Base34428429 _ssavalue_0;
     _ssavalue_0 = x428;
     return map_226(_ssavalue_0, ts);
 }
 // NTuple{6,Tuple{Float32}}
 struct  __attribute__ ((packed)) TYPTuple_float_float_float_float_float_float{
     float field1;
     float field2;
     float field3;
     float field4;
     float field5;
     float field6;
 };
 typedef struct TYPTuple_float_float_float_float_float_float Tuple_float_float_float_float_float_float;

 // NTuple{5,Tuple{Float32}}
 struct  __attribute__ ((packed)) TYPTuple_float_float_float_float_float{
     float field1;
     float field2;
     float field3;
     float field4;
     float field5;
 };
 typedef struct TYPTuple_float_float_float_float_float Tuple_float_float_float_float_float;

 // NTuple{4,Tuple{Float32}}
 struct  __attribute__ ((packed)) TYPTuple_float_float_float_float{
     float field1;
     float field2;
     float field3;
     float field4;
 };
 typedef struct TYPTuple_float_float_float_float Tuple_float_float_float_float;

 // Tuple{Tuple{Float32},Tuple{Float32},Tuple{Float32}}
 struct  __attribute__ ((packed)) TYPTuple_float_float_float{
     float field1;
     float field2;
     float field3;
 };
 typedef struct TYPTuple_float_float_float Tuple_float_float_float;

 // (map, Tuple{Base.#tail,Tuple{Tuple{Float32,Float32},Tuple{Float32,Float32},Tuple{Float32,Float32}}})
 Tuple_float_float_float map_227(Base34tail f, Tuple_float2_float2_float2 t)
 {
     return (Tuple_float_float_float){tail_57(t.field1), tail_57(t.field2), tail_57(t.field3)};
 }
 // (map, Tuple{Base.#tail,NTuple{4,Tuple{Float32,Float32}}})
 Tuple_float_float_float_float map_228(Base34tail f, Tuple_float2_float2_float2_float2 t)
 {
     Tuple_float_float_float x44_apply_tmp41581;
     x44_apply_tmp41581 = map_227(f, tail_219(t));
     return (Tuple_float_float_float_float){tail_57(t.field1), x44_apply_tmp41581.field1, x44_apply_tmp41581.field2, x44_apply_tmp41581.field3};
 }
 // (map, Tuple{Base.#tail,NTuple{5,Tuple{Float32,Float32}}})
 Tuple_float_float_float_float_float map_229(Base34tail f, Tuple_float2_float2_float2_float2_float2 t)
 {
     Tuple_float_float_float_float x44_apply_tmp41580;
     x44_apply_tmp41580 = map_228(f, tail_222(t));
     return (Tuple_float_float_float_float_float){tail_57(t.field1), x44_apply_tmp41580.field1, x44_apply_tmp41580.field2, x44_apply_tmp41580.field3, x44_apply_tmp41580.field4};
 }
 // (map, Tuple{Base.#tail,NTuple{6,Tuple{Float32,Float32}}})
 Tuple_float_float_float_float_float_float map_230(Base34tail f, Tuple_float2_float2_float2_float2_float2_float2 t)
 {
     Tuple_float_float_float_float_float x44_apply_tmp41579;
     x44_apply_tmp41579 = map_229(f, tail_225(t));
     return (Tuple_float_float_float_float_float_float){tail_57(t.field1), x44_apply_tmp41579.field1, x44_apply_tmp41579.field2, x44_apply_tmp41579.field3, x44_apply_tmp41579.field4, x44_apply_tmp41579.field5};
 }
 // (Base.tails, NTuple{6,Tuple{Float32,Float32}})
 Tuple_float_float_float_float_float_float tails_224(Tuple_float2_float2_float2_float2_float2_float2 ts)
 {
     return map_230(FUNC_INST_Base34tail, ts);
 }
 // (map, Tuple{Base.##28#29,Tuple{Tuple{Float32},Tuple{Float32},Tuple{Float32}}})
 float3 map_231(Base34428429 f, Tuple_float_float_float t)
 {
     return (float3){x428_107(t.field1), x428_107(t.field2), x428_107(t.field3)};
 }
 // (Base.argtail, NTuple{4,Tuple{Float32}})
 Tuple_float_float_float argtail_232(float x, Tuple_float_float_float rest)
 {
     return rest;
 }
 // (Base.tail, Tuple{NTuple{4,Tuple{Float32}}})
 Tuple_float_float_float tail_233(Tuple_float_float_float_float x)
 {
     Tuple_float_float_float_float x44_apply_tmp41590;
     x44_apply_tmp41590 = x;
     return argtail_232(x44_apply_tmp41590.field1, (Tuple_float_float_float){x44_apply_tmp41590.field2, x44_apply_tmp41590.field3, x44_apply_tmp41590.field4});
 }
 // (map, Tuple{Base.##28#29,NTuple{4,Tuple{Float32}}})
 float4 map_234(Base34428429 f, Tuple_float_float_float_float t)
 {
     float3 x44_apply_tmp41589;
     x44_apply_tmp41589 = map_231(f, tail_233(t));
     return (float4){x428_107(t.field1), x44_apply_tmp41589.s0, x44_apply_tmp41589.s1, x44_apply_tmp41589.s2};
 }
 // (Base.argtail, NTuple{5,Tuple{Float32}})
 Tuple_float_float_float_float argtail_235(float x, Tuple_float_float_float_float rest)
 {
     return rest;
 }
 // (Base.tail, Tuple{NTuple{5,Tuple{Float32}}})
 Tuple_float_float_float_float tail_236(Tuple_float_float_float_float_float x)
 {
     Tuple_float_float_float_float_float x44_apply_tmp41591;
     x44_apply_tmp41591 = x;
     return argtail_235(x44_apply_tmp41591.field1, (Tuple_float_float_float_float){x44_apply_tmp41591.field2, x44_apply_tmp41591.field3, x44_apply_tmp41591.field4, x44_apply_tmp41591.field5});
 }
 // (map, Tuple{Base.##28#29,NTuple{5,Tuple{Float32}}})
 Tuple_float_float_float_float_float map_237(Base34428429 f, Tuple_float_float_float_float_float t)
 {
     float4 x44_apply_tmp41588;
     x44_apply_tmp41588 = map_234(f, tail_236(t));
     return (Tuple_float_float_float_float_float){x428_107(t.field1), x44_apply_tmp41588.s0, x44_apply_tmp41588.s1, x44_apply_tmp41588.s2, x44_apply_tmp41588.s3};
 }
 // (Base.argtail, NTuple{6,Tuple{Float32}})
 Tuple_float_float_float_float_float argtail_238(float x, Tuple_float_float_float_float_float rest)
 {
     return rest;
 }
 // (Base.tail, Tuple{NTuple{6,Tuple{Float32}}})
 Tuple_float_float_float_float_float tail_239(Tuple_float_float_float_float_float_float x)
 {
     Tuple_float_float_float_float_float_float x44_apply_tmp41592;
     x44_apply_tmp41592 = x;
     return argtail_238(x44_apply_tmp41592.field1, (Tuple_float_float_float_float_float){x44_apply_tmp41592.field2, x44_apply_tmp41592.field3, x44_apply_tmp41592.field4, x44_apply_tmp41592.field5, x44_apply_tmp41592.field6});
 }
 // (map, Tuple{Base.##28#29,NTuple{6,Tuple{Float32}}})
 Tuple_float_float_float_float_float_float map_240(Base34428429 f, Tuple_float_float_float_float_float_float t)
 {
     Tuple_float_float_float_float_float x44_apply_tmp41587;
     x44_apply_tmp41587 = map_237(f, tail_239(t));
     return (Tuple_float_float_float_float_float_float){x428_107(t.field1), x44_apply_tmp41587.field1, x44_apply_tmp41587.field2, x44_apply_tmp41587.field3, x44_apply_tmp41587.field4, x44_apply_tmp41587.field5};
 }
 // (Base.heads, NTuple{6,Tuple{Float32}})
 Tuple_float_float_float_float_float_float heads_238(Tuple_float_float_float_float_float_float ts)
 {
     Base34428429 x428;
     x428 = (Base34428429){0.0f};
     Base34428429 _ssavalue_0;
     _ssavalue_0 = x428;
     return map_240(_ssavalue_0, ts);
 }
 // NTuple{6,Tuple{}}
 struct  __attribute__ ((packed)) TYPTuple_EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple_{
     EmptyTuple_ field1;
     EmptyTuple_ field2;
     EmptyTuple_ field3;
     EmptyTuple_ field4;
     EmptyTuple_ field5;
     EmptyTuple_ field6;
 };
 typedef struct TYPTuple_EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple_ Tuple_EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple_;

 // NTuple{5,Tuple{}}
 struct  __attribute__ ((packed)) TYPTuple_EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple_{
     EmptyTuple_ field1;
     EmptyTuple_ field2;
     EmptyTuple_ field3;
     EmptyTuple_ field4;
     EmptyTuple_ field5;
 };
 typedef struct TYPTuple_EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple_ Tuple_EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple_;

 // NTuple{4,Tuple{}}
 struct  __attribute__ ((packed)) TYPTuple_EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple_{
     EmptyTuple_ field1;
     EmptyTuple_ field2;
     EmptyTuple_ field3;
     EmptyTuple_ field4;
 };
 typedef struct TYPTuple_EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple_ Tuple_EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple_;

 // (map, Tuple{Base.#tail,Tuple{Tuple{Float32},Tuple{Float32},Tuple{Float32}}})
 Tuple_EmptyTuple__EmptyTuple__EmptyTuple_ map_241(Base34tail f, Tuple_float_float_float t)
 {
     return (Tuple_EmptyTuple__EmptyTuple__EmptyTuple_){tail_107(t.field1), tail_107(t.field2), tail_107(t.field3)};
 }
 // (map, Tuple{Base.#tail,NTuple{4,Tuple{Float32}}})
 Tuple_EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple_ map_242(Base34tail f, Tuple_float_float_float_float t)
 {
     Tuple_EmptyTuple__EmptyTuple__EmptyTuple_ x44_apply_tmp41595;
     x44_apply_tmp41595 = map_241(f, tail_233(t));
     return (Tuple_EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple_){tail_107(t.field1), x44_apply_tmp41595.field1, x44_apply_tmp41595.field2, x44_apply_tmp41595.field3};
 }
 // (map, Tuple{Base.#tail,NTuple{5,Tuple{Float32}}})
 Tuple_EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple_ map_243(Base34tail f, Tuple_float_float_float_float_float t)
 {
     Tuple_EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple_ x44_apply_tmp41594;
     x44_apply_tmp41594 = map_242(f, tail_236(t));
     return (Tuple_EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple_){tail_107(t.field1), x44_apply_tmp41594.field1, x44_apply_tmp41594.field2, x44_apply_tmp41594.field3, x44_apply_tmp41594.field4};
 }
 // (map, Tuple{Base.#tail,NTuple{6,Tuple{Float32}}})
 Tuple_EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple_ map_244(Base34tail f, Tuple_float_float_float_float_float_float t)
 {
     Tuple_EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple_ x44_apply_tmp41593;
     x44_apply_tmp41593 = map_243(f, tail_239(t));
     return (Tuple_EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple_){tail_107(t.field1), x44_apply_tmp41593.field1, x44_apply_tmp41593.field2, x44_apply_tmp41593.field3, x44_apply_tmp41593.field4, x44_apply_tmp41593.field5};
 }
 // (Base.tails, NTuple{6,Tuple{Float32}})
 Tuple_EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple_ tails_238(Tuple_float_float_float_float_float_float ts)
 {
     return map_244(FUNC_INST_Base34tail, ts);
 }
 // (map, Tuple{##22#26,Tuple{},Tuple{},Tuple{},Tuple{},Tuple{},Tuple{}})
 EmptyTuple_ map_245(x4422426 f, Tuple_EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple__EmptyTuple_ x4unused4)
 {
     return (EmptyTuple_){0.0f};
 }
 // (map, Tuple{##22#26,Tuple{Float32},Tuple{Float32},Tuple{Float32},Tuple{Float32},Tuple{Float32},Tuple{Float32}})
 float map_246(x4422426 f, float t1, float t2, Tuple_float_float_float_float ts)
 {
     Tuple_float_float_float_float_float_float x44_apply_tmp41583;
     x44_apply_tmp41583 = Sugar.InlineNode(Any[:(_6::NTuple{4,Tuple{Float32}}::NTuple{4,Tuple{Float32}}), :(_6::NTuple{4,Tuple{Float32}} = _5::NTuple{4,Tuple{Float32}})], :((Base.heads)((NTuple{6,Tuple{Float32}}){_3::Tuple{Float32}, _4::Tuple{Float32}, (getfield)(_6::NTuple{4,Tuple{Float32}}, field1)::Tuple{Float32}, (getfield)(_6::NTuple{4,Tuple{Float32}}, field2)::Tuple{Float32}, (getfield)(_6::NTuple{4,Tuple{Float32}}, field3)::Tuple{Float32}, (getfield)(_6::NTuple{4,Tuple{Float32}}, field4)::Tuple{Float32}}::NTuple{6,Tuple{Float32}})::NTuple{6,Float32}));
     EmptyTuple_ x44_apply_tmp41586;
     x44_apply_tmp41586 = Sugar.InlineNode(Any[:(_9::NTuple{6,Tuple{}}::NTuple{6,Tuple{}}), :(_9::NTuple{6,Tuple{}} = Sugar.InlineNode(Any[:(_8::NTuple{4,Tuple{Float32}}::NTuple{4,Tuple{Float32}}), :(_8::NTuple{4,Tuple{Float32}} = _5::NTuple{4,Tuple{Float32}})], :((Base.tails)((NTuple{6,Tuple{Float32}}){_3::Tuple{Float32}, _4::Tuple{Float32}, (getfield)(_8::NTuple{4,Tuple{Float32}}, field1)::Tuple{Float32}, (getfield)(_8::NTuple{4,Tuple{Float32}}, field2)::Tuple{Float32}, (getfield)(_8::NTuple{4,Tuple{Float32}}, field3)::Tuple{Float32}, (getfield)(_8::NTuple{4,Tuple{Float32}}, field4)::Tuple{Float32}}::NTuple{6,Tuple{Float32}})::NTuple{6,Tuple{}})))], :((map)(_2::##22#26, (NTuple{6,Tuple{}}){(getfield)(_9::NTuple{6,Tuple{}}, field1)::Tuple{}, (getfield)(_9::NTuple{6,Tuple{}}, field2)::Tuple{}, (getfield)(_9::NTuple{6,Tuple{}}, field3)::Tuple{}, (getfield)(_9::NTuple{6,Tuple{}}, field4)::Tuple{}, (getfield)(_9::NTuple{6,Tuple{}}, field5)::Tuple{}, (getfield)(_9::NTuple{6,Tuple{}}, field6)::Tuple{}}::NTuple{6,Tuple{}})::Tuple{}));
     return (float){x422_212(x44_apply_tmp41583.field1, x44_apply_tmp41583.field2, x44_apply_tmp41583.field3, x44_apply_tmp41583.field4, x44_apply_tmp41583.field5, x44_apply_tmp41583.field6)};
 }
 // (map, Tuple{##22#26,Tuple{Float32,Float32},Tuple{Float32,Float32},Tuple{Float32,Float32},Tuple{Float32,Float32},Tuple{Float32,Float32},Tuple{Float32,Float32}})
 float2 map_247(x4422426 f, float2 t1, float2 t2, Tuple_float2_float2_float2_float2 ts)
 {
     Tuple_float_float_float_float_float_float x44_apply_tmp41569;
     x44_apply_tmp41569 = Sugar.InlineNode(Any[:(_6::NTuple{4,Tuple{Float32,Float32}}::NTuple{4,Tuple{Float32,Float32}}), :(_6::NTuple{4,Tuple{Float32,Float32}} = _5::NTuple{4,Tuple{Float32,Float32}})], :((Base.heads)((NTuple{6,Tuple{Float32,Float32}}){_3::Tuple{Float32,Float32}, _4::Tuple{Float32,Float32}, (getfield)(_6::NTuple{4,Tuple{Float32,Float32}}, field1)::Tuple{Float32,Float32}, (getfield)(_6::NTuple{4,Tuple{Float32,Float32}}, field2)::Tuple{Float32,Float32}, (getfield)(_6::NTuple{4,Tuple{Float32,Float32}}, field3)::Tuple{Float32,Float32}, (getfield)(_6::NTuple{4,Tuple{Float32,Float32}}, field4)::Tuple{Float32,Float32}}::NTuple{6,Tuple{Float32,Float32}})::NTuple{6,Float32}));
     float x44_apply_tmp41572;
     x44_apply_tmp41572 = Sugar.InlineNode(Any[:(_9::NTuple{6,Tuple{Float32}}::NTuple{6,Tuple{Float32}}), :(_9::NTuple{6,Tuple{Float32}} = Sugar.InlineNode(Any[:(_8::NTuple{4,Tuple{Float32,Float32}}::NTuple{4,Tuple{Float32,Float32}}), :(_8::NTuple{4,Tuple{Float32,Float32}} = _5::NTuple{4,Tuple{Float32,Float32}})], :((Base.tails)((NTuple{6,Tuple{Float32,Float32}}){_3::Tuple{Float32,Float32}, _4::Tuple{Float32,Float32}, (getfield)(_8::NTuple{4,Tuple{Float32,Float32}}, field1)::Tuple{Float32,Float32}, (getfield)(_8::NTuple{4,Tuple{Float32,Float32}}, field2)::Tuple{Float32,Float32}, (getfield)(_8::NTuple{4,Tuple{Float32,Float32}}, field3)::Tuple{Float32,Float32}, (getfield)(_8::NTuple{4,Tuple{Float32,Float32}}, field4)::Tuple{Float32,Float32}}::NTuple{6,Tuple{Float32,Float32}})::NTuple{6,Tuple{Float32}})))], :((map)(_2::##22#26, (getfield)(_9::NTuple{6,Tuple{Float32}}, field1)::Tuple{Float32}, (getfield)(_9::NTuple{6,Tuple{Float32}}, field2)::Tuple{Float32}, (NTuple{4,Tuple{Float32}}){(getfield)(_9::NTuple{6,Tuple{Float32}}, field3)::Tuple{Float32}, (getfield)(_9::NTuple{6,Tuple{Float32}}, field4)::Tuple{Float32}, (getfield)(_9::NTuple{6,Tuple{Float32}}, field5)::Tuple{Float32}, (getfield)(_9::NTuple{6,Tuple{Float32}}, field6)::Tuple{Float32}}::NTuple{4,Tuple{Float32}})::Tuple{Float32}));
     return (float2){x422_212(x44_apply_tmp41569.field1, x44_apply_tmp41569.field2, x44_apply_tmp41569.field3, x44_apply_tmp41569.field4, x44_apply_tmp41569.field5, x44_apply_tmp41569.field6), x44_apply_tmp41572};
 }
 // (map, Tuple{##22#26,Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32}})
 float3 map_248(x4422426 f, float3 t1, float3 t2, Tuple_float3_float3_float3_float3 ts)
 {
     Tuple_float_float_float_float_float_float x44_apply_tmp41555;
     x44_apply_tmp41555 = Sugar.InlineNode(Any[:(_6::NTuple{4,Tuple{Float32,Float32,Float32}}::NTuple{4,Tuple{Float32,Float32,Float32}}), :(_6::NTuple{4,Tuple{Float32,Float32,Float32}} = _5::NTuple{4,Tuple{Float32,Float32,Float32}})], :((Base.heads)((NTuple{6,Tuple{Float32,Float32,Float32}}){_3::Tuple{Float32,Float32,Float32}, _4::Tuple{Float32,Float32,Float32}, (getfield)(_6::NTuple{4,Tuple{Float32,Float32,Float32}}, field1)::Tuple{Float32,Float32,Float32}, (getfield)(_6::NTuple{4,Tuple{Float32,Float32,Float32}}, field2)::Tuple{Float32,Float32,Float32}, (getfield)(_6::NTuple{4,Tuple{Float32,Float32,Float32}}, field3)::Tuple{Float32,Float32,Float32}, (getfield)(_6::NTuple{4,Tuple{Float32,Float32,Float32}}, field4)::Tuple{Float32,Float32,Float32}}::NTuple{6,Tuple{Float32,Float32,Float32}})::NTuple{6,Float32}));
     float2 x44_apply_tmp41558;
     x44_apply_tmp41558 = Sugar.InlineNode(Any[:(_9::NTuple{6,Tuple{Float32,Float32}}::NTuple{6,Tuple{Float32,Float32}}), :(_9::NTuple{6,Tuple{Float32,Float32}} = Sugar.InlineNode(Any[:(_8::NTuple{4,Tuple{Float32,Float32,Float32}}::NTuple{4,Tuple{Float32,Float32,Float32}}), :(_8::NTuple{4,Tuple{Float32,Float32,Float32}} = _5::NTuple{4,Tuple{Float32,Float32,Float32}})], :((Base.tails)((NTuple{6,Tuple{Float32,Float32,Float32}}){_3::Tuple{Float32,Float32,Float32}, _4::Tuple{Float32,Float32,Float32}, (getfield)(_8::NTuple{4,Tuple{Float32,Float32,Float32}}, field1)::Tuple{Float32,Float32,Float32}, (getfield)(_8::NTuple{4,Tuple{Float32,Float32,Float32}}, field2)::Tuple{Float32,Float32,Float32}, (getfield)(_8::NTuple{4,Tuple{Float32,Float32,Float32}}, field3)::Tuple{Float32,Float32,Float32}, (getfield)(_8::NTuple{4,Tuple{Float32,Float32,Float32}}, field4)::Tuple{Float32,Float32,Float32}}::NTuple{6,Tuple{Float32,Float32,Float32}})::NTuple{6,Tuple{Float32,Float32}})))], :((map)(_2::##22#26, (getfield)(_9::NTuple{6,Tuple{Float32,Float32}}, field1)::Tuple{Float32,Float32}, (getfield)(_9::NTuple{6,Tuple{Float32,Float32}}, field2)::Tuple{Float32,Float32}, (NTuple{4,Tuple{Float32,Float32}}){(getfield)(_9::NTuple{6,Tuple{Float32,Float32}}, field3)::Tuple{Float32,Float32}, (getfield)(_9::NTuple{6,Tuple{Float32,Float32}}, field4)::Tuple{Float32,Float32}, (getfield)(_9::NTuple{6,Tuple{Float32,Float32}}, field5)::Tuple{Float32,Float32}, (getfield)(_9::NTuple{6,Tuple{Float32,Float32}}, field6)::Tuple{Float32,Float32}}::NTuple{4,Tuple{Float32,Float32}})::Tuple{Float32,Float32}));
     return (float3){x422_212(x44_apply_tmp41555.field1, x44_apply_tmp41555.field2, x44_apply_tmp41555.field3, x44_apply_tmp41555.field4, x44_apply_tmp41555.field5, x44_apply_tmp41555.field6), x44_apply_tmp41558.s0, x44_apply_tmp41558.s1};
 }
 // (broadcast, Tuple{##22#26,Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32}})
 float3 broadcast_248(x4422426 f, float3 t, Tuple_float3_float3_float3_float3_float3 ts)
 {
     Tuple_float3_float3_float3_float3_float3 x44_apply_tmp41553;
     x44_apply_tmp41553 = ts;
     return map_248(f, t, x44_apply_tmp41553.field1, (Tuple_float3_float3_float3_float3){x44_apply_tmp41553.field2, x44_apply_tmp41553.field3, x44_apply_tmp41553.field4, x44_apply_tmp41553.field5});
 }
 // (staggered_velocity, Tuple{CLArrays.DeviceArray{Tuple{Float32,Float32,Float32},3,Transpiler.CLIntrinsics.GlobalPointer{Tuple{Float32,Float32,Float32}}},Tuple{Float32,Float32,Float32},Tuple{Float32,Float32,Float32},Tuple{UInt32,UInt32,UInt32},Tuple{UInt32,UInt32,UInt32}})
 float3 staggered_velocity_249(DeviceArray_float3_3___global1float3121 velocity, float3 point, float3 d, uint3 gs, uint3 res)
 {
     float3 w2;
     float3 w1;
     float3 w;
     float3 pp2;
     float3 pp;
     float3 vn;
     float3 pzp;
     float3 pyp;
     float3 pxp;
     float3 v0;
     uint3 ip;
     uint3 i;
     float3 p;
     x4422426 x422;
     x4421425 x421;
     x4420424 x420;
     x4419423 x419;
     p = broadcast_157(FUNC_INST_Base34mod, point, (Tuple_uint3){gs});
     x419 = (x4419423){0.0f};
     x4419423 _ssavalue_0;
     _ssavalue_0 = x419;
     float3 _ssavalue_1;
     _ssavalue_1 = p;
     float3 _ssavalue_2;
     _ssavalue_2 = d;
     uint _ssavalue_3;
     _ssavalue_3 = (uint)(1);
     i = broadcast_175(_ssavalue_0, _ssavalue_1, (Tuple_float3_uint){_ssavalue_2, _ssavalue_3});
     x420 = (x4420424){0.0f};
     x4420424 _ssavalue_4;
     _ssavalue_4 = x420;
     uint3 _ssavalue_5;
     _ssavalue_5 = i;
     uint3 _ssavalue_6;
     _ssavalue_6 = res;
     uint _ssavalue_7;
     _ssavalue_7 = (uint)(1);
     ip = broadcast_187(_ssavalue_4, _ssavalue_5, (Tuple_uint3_uint){_ssavalue_6, _ssavalue_7});
     v0 = getindex_135(velocity, (uint3){i.s0, i.s1, i.s2});
     pxp = getindex_135(velocity, (uint3){ip.s0, i.s1, i.s2});
     pyp = getindex_135(velocity, (uint3){i.s0, ip.s1, i.s2});
     pzp = getindex_135(velocity, (uint3){i.s0, i.s1, ip.s2});
     vn = (float3){getindex_135(velocity, (uint3){i.s0, ip.s1, ip.s2}).s0, getindex_135(velocity, (uint3){ip.s0, i.s1, ip.s2}).s1, getindex_135(velocity, (uint3){ip.s0, ip.s1, i.s2}).s2};
     pp = (float3){pyp.s0, pxp.s1, pxp.s2};
     pp2 = (float3){pzp.s0, pzp.s1, pyp.s2};
     x421 = (x4421425){0.0f};
     x4421425 _ssavalue_8;
     _ssavalue_8 = x421;
     float3 _ssavalue_9;
     _ssavalue_9 = p;
     uint3 _ssavalue_10;
     _ssavalue_10 = i;
     w = broadcast_201(_ssavalue_8, _ssavalue_9, (Tuple_uint3_float3){_ssavalue_10, d});
     w1 = (float3){w.s2, w.s2, w.s1};
     w2 = (float3){w.s1, w.s0, w.s0};
     x422 = (x4422426){0.0f};
     x4422426 _ssavalue_11;
     _ssavalue_11 = x422;
     float3 _ssavalue_12;
     _ssavalue_12 = w1;
     float3 _ssavalue_13;
     _ssavalue_13 = w2;
     float3 _ssavalue_14;
     _ssavalue_14 = v0;
     float3 _ssavalue_15;
     _ssavalue_15 = pp;
     float3 _ssavalue_16;
     _ssavalue_16 = pp2;
     return broadcast_248(_ssavalue_11, _ssavalue_12, (Tuple_float3_float3_float3_float3_float3){_ssavalue_13, _ssavalue_14, _ssavalue_15, _ssavalue_16, vn});
 }
 // CLArrays.KernelState
 struct  __attribute__ ((packed)) TYPKernelState{
     int empty;
 };
 typedef struct TYPKernelState KernelState;

 // ##51#52
 __constant int FUNC_INST_x4451452 = 0;
 typedef int x4451452; // empty type emitted as an int
 // ########################
 // Main inner function
 // (#51, (CLArrays.KernelState, CLArrays.DeviceArray{Tuple{Float32,Float32,Float32},3,Transpiler.CLIntrinsics.GlobalPointer{Tuple{Float32,Float32,Float32}}}, Float32, Tuple{Float32,Float32,Float32}, Tuple{UInt32,UInt32,UInt32}, Tuple{UInt32,UInt32,UInt32}))
 __kernel float3 x451_250(KernelState state, DeviceArray_float3_3_HostPtr_float3 x44velo41500, float dt, float3 d, uint3 ps, uint3 gr, __global float3 *  x44ptr_1_velo141501){
     float3 point;
     DeviceArray_float3_3___global1float3121 velo;
     velo = reconstruct_34(x44velo41500, x44ptr_1_velo141501);
     point = (float3){0.5f, 0.5f, 0.5f};
     return staggered_velocity_249(velo, point, d, ps, gr);
 }
With following build error:
<kernel>:842:34: error: redefinition of 'TYPTuple_float_uint_float'
struct  __attribute__ ((packed)) TYPTuple_float_uint_float{
                                 ^
<kernel>:732:34: note: previous definition is here
struct  __attribute__ ((packed)) TYPTuple_float_uint_float{
                                 ^
