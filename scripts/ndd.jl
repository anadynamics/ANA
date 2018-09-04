#!/usr/bin/env julia
###############################################################################
# Utility to displace a PDB along many vectors, generating 1 PDB for each
# vector, useful for NDD analysis.
# Calpha only.
# by https://github.com/pgbarletta
###############################################################################
using Chemfiles
using StaticArrays
using ArgParse
##########
# functions
##########
function read_ptraj_modes(filename, nmodes::Int64=0, norma::Bool=true)
    modes_text = readdlm(filename, skipstart=0, skipblanks=true, comments=true,
        comment_char='\*')

    if nmodes == 0
        nmodes = modes_text[1, 5]
    end
    modes_elements = modes_text[2, 1]

    ncoords = convert(Int64, modes_elements)
    lines = ceil(Int64, ncoords/7)
    rest = convert(Int64, ncoords % 7)

    eval = Array{Float64}(nmodes);
    mode = Array{Float64}(ncoords, nmodes);
    temp1 = Array{Float64}(ncoords, 1);
    temp2 = Array{Float64}(ncoords+(7-rest));

    j=lines + 1 + 2 # 1 p/ q lea la prox linea 2 por el header
    for i=1:nmodes
        eval[i] = modes_text[j, 2]
        temp = permutedims(modes_text[(j+1):(lines+j), :], [2, 1])
        temp2 = reshape(temp, ncoords+(7-rest))
        for k=(rest+1):7
            pop!(temp2)
        end
        mode[:, i] = temp2
        j = j + lines + 1
    end

    if norma == true
        for i=1:nmodes
            mode[: ,i] = mode[:, i] / norm(mode[:, i])
        end
    end

    return mode, eval
end
#########
function displaceAA(in_frm, aa, aa_3, in_vec)
    # Preparo variables
    const in_top = Topology(in_frm)
    natoms = convert(Int64, size(in_top))
    const in_xyz = positions(in_frm)

    # Determino orden de residuos (hay q actualizar el Julia Chemfiles)
    tmp = Array{Int64}(aa)
    ids = Array{Int64}(aa)
    [ ids[i+1] = convert(Int64, id((Residue(in_top, i)))) for i = 0:aa-1 ]
    const idx = sortperm(ids)
    # Determino el nro de atomos de c/ aminoácido
    [ tmp[i+1] = size(Residue(in_top, i)) for i = 0:aa-1 ]
    const natom_aa = tmp[idx]

    # Paso el vector columna de tamaño 1xaa_3 a 3xaa
    const vector = reshape(in_vec, 3, aa)
    # Adapto el vector p/ darle la misma forma q la matriz de coordenadas
    sum_mat = Array{Float64}(3, natoms)
    cursor = 0
    for i = 1:aa
        if i == 1
            sum_mat[:, 1:natom_aa[i]] = repmat(vector[:, 1], 1, natom_aa[i])
            cursor = natom_aa[i]
            continue
        end
        const rango = collect(cursor+1:cursor + natom_aa[i])
        sum_mat[:, rango] = repmat(vector[:, i], 1, natom_aa[i])
        cursor += natom_aa[i]
    end

    # Listo, ahora puedo mover el pdb
    out_frm = deepcopy(in_frm)
    out_xyz = positions(out_frm)

    # Tengo q hacer esto por ahora. Hasta q arreglemos Chemfiles.
    for i = 1:size(in_xyz)[1]
        for j = 1:size(in_xyz)[2]
            out_xyz[i, j]  = round(in_xyz[i, j] + sum_mat[i, j], 3)
        end
    end
    return out_frm
end
#########
function displaceAtoms(mod_pdb, vector1, multiplier)
  # Preparo variables
    pdb = copy(mod_pdb)
    struct_xyz = coordinatesmatrix(pdb)
#    new_struct_xyz = copy(struct_xyz)
    vector = Array{Float64}(1, 3)

    # Adapto el vector p/ darle la misma forma q la matriz de coordenadas
    for i = 1:3:length(vector1)
        if i== 1
            vector = reshape(vector1[i:i+2], 1, 3)
            continue
        end
        vector = vcat(vector, reshape(vector1[i:i+2], 1, 3))
    end

    # Listo, ahora puedo mover el pdb
    new_struct_xyz  = struct_xyz + vector .* multiplier
    pdb = change_coordinates(pdb, new_struct_xyz);
   return pdb
end
#########
# Arg Parse settings
s = ArgParseSettings()
@add_arg_table s begin
    "--in_pdb_filename", "-p"
        help = "Input PDB."
        arg_type = String
        required = true
    "--modes_filename", "-v"
        help = "Input modes."
        arg_type = String
        required = true
    "--mul", "-m"
        help = "Multiplier."
        arg_type = Int
        required = true
    "--suffix", "-o"
        help = "Output PDBs suffix"
        arg_type = String
        required = true
    "--amber_modes", "-a"
        help = "Mark true when reading from Amber PCA modes kind of file. Default: false."
        arg_type = Bool
        required = false
        default = false
    "--weights_filename", "-w"
        help = "Input 1/weights, if desired. Default: none"
        arg_type = String
        required = false
        default = "none"
end

##########
# main program
##########

# Read arguments from console
parsed_args = parse_args(ARGS, s)
args = Array{Any, 1}(0)
for (arg, val) in parsed_args
    arg = Symbol(arg)
    @eval (($arg) = ($val))
end
println("Input parameters:")
println("INPDB          ", in_pdb_filename)
println("MODES          ", modes_filename)
println("MUL            ", mul)
println("SUFFIX         ", suffix)
println("AMBER_MODES    ", amber_modes)
println("WEIGHTS        ", weights_filename)

# Read PDB
const in_trj = Trajectory(in_pdb_filename)
const in_frm = read(in_trj)
const in_top = Topology(in_frm)
const aa = convert(Int64, count_residues(in_top))
const aa_3 = aa * 3
# Initialize modes vectors
in_modes = Array{Float64, 2}(aa_3, aa)
in_evals = Array{Float64, 1}(aa_3)
if (amber_modes)
    try
        in_modes, in_evals = read_ptraj_modes(modes_filename)
    catch e
        println(modes_filename, " error:")
        error(e)
    end
else
    try
        in_modes = convert(Array{Float64, 2}, readdlm(modes_filename))
    catch e
        println(modes_filename, " error:")
        error(e)
    end
end
const l_modes = size(in_modes)[1]
const n_modes = size(in_modes)[2]
# Check input modes
if l_modes != aa_3
# Asumo q el modo es de Calpha y está ordenado en una columna
    error(string("\n\n", " Input PDB has ", aa, " residues, so input ",
        "modes should be ", aa_3, " long."))
end

# Initialize weights
weights = MVector{n_modes, Float64}()
if (amber_modes)
    if weights_filename != "none"
        try
            weights = convert(MVector{n_modes, Float64},
                readdlm(weights_filename)[:, 1])
        catch e
            println(weights_filename, " error:")
            error(e)
        end
    else
        try
            weights = convert(MVector{n_modes, Float64}, in_evals)
        catch e
            println("Number of amber PCA eigenvalues does not match number of ",
                "PCA modes")
            error(e)
        end
    end
else
    if weights_filename != "none"
        try
            weights = convert(MVector{n_modes, Float64},
                readdlm(weights_filename)[:, 1])
        catch e
            println(weights_filename, " error:")
            error(e)
        end
    else
        weights = convert(MVector{n_modes, Float64}, fill(1.0, n_modes))
    end
end

# Ahora desplazo
pdb_names = Array{String}(n_modes)
for i = 1:n_modes
    # Escalo vector
    const modo = in_modes[:, i] .* mul ./ weights[i]
    # Desplazo
    const out_frm = displaceAA(in_frm, aa, aa_3, modo);
    # Y guardo
    pdb_names[i] = joinpath(pwd(), string(i, "_", suffix, ".pdb"))
    out_trj = Trajectory(pdb_names[i], 'w')
    write(out_trj, out_frm)
end
# Write in_ndd file
const in_ndd_filename = string("in_ndd_",
    splitext(basename(in_pdb_filename))[1])
const in_ndd_location = dirname(in_pdb_filename)
writedlm(joinpath(in_ndd_location, in_ndd_filename), pdb_names)
