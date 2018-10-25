#!/usr/bin/env julia
###############################################################################
# Utility to displace a PDB along many vectors, generating 1 PDB for each
# vector, useful for NDD analysis.
# Calpha only.
# by https://github.com/pgbarletta
###############################################################################
using Chemfiles, JUMD
using StaticArrays
using ArgParse
using DelimitedFiles, LinearAlgebra, FileIO

function displaceAA(in_frm, aa, aa_3, in_vec)
    # Preparo variables
    in_top = Topology(in_frm)
    natoms = convert(Int64, size(in_top))
    in_xyz = positions(in_frm)

    # Determino orden de residuos (hay q actualizar el Julia Chemfiles)
    tmp = Array{Int64}(undef, aa)
    ids = Array{Int64}(undef, aa)
    [ ids[i+1] = convert(Int64, id((Residue(in_top, i)))) for i = 0:aa-1 ]
    idx = sortperm(ids)
    # Determino el nro de atomos de c/ aminoácido
    [ tmp[i+1] = size(Residue(in_top, i)) for i = 0:aa-1 ]
    natom_aa = tmp[idx]

    # Paso el vector columna de tamaño 1xaa_3 a 3xaa
    vector = reshape(in_vec, 3, aa)
    # Adapto el vector p/ darle la misma forma q la matriz de coordenadas
    sum_mat = Array{Float64}(undef, 3, natoms)
    cursor = 0
    for i = 1:aa
        if i == 1
            sum_mat[:, 1:natom_aa[i]] = repmat(vector[:, 1], 1, natom_aa[i])
            cursor = natom_aa[i]
            continue
        end
        rango = collect(cursor+1:cursor + natom_aa[i])
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
args = Array{Any, 1}(undef, 0)
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
in_modes = Array{Float64, 2}(undef, aa_3, aa)
in_evals = Array{Float64, 1}(undef, aa_3)
if (amber_modes)
    try
        in_modes, in_evals = JUMD.readPtrajModes(modes_filename)
    catch e
        println("")
        println("Error when reading Amber modes: ", modes_filename)
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
pdb_names = Array{String}(undef, n_modes)
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
