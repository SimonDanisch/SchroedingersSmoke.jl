packages = (
    "GLVisualize",
    "GLAbstraction",
    "GLWindow",
    "GeometryTypes",
    "MeshIO",
)

for pkg in packages
    Pkg.checkout(pkg, "sd/staticarrays")
end
