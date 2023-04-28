# module NasReader

# num_nodes_per_cell = Dict(
#     "vertex"=> 1,
#     "line"=> 2,
#     "triangle"=> 3,
#     "quad"=> 4,
#     "quad8"=> 8,
#     "tetra"=> 4,
#     "hexahedron"=> 8,
#     "hexahedron20"=> 20,
#     "hexahedron24"=> 24,
#     "wedge"=> 6,
#     "pyramid"=> 5,
#     #
#     "line3"=> 3,
#     "triangle6"=> 6,
#     "quad9"=> 9,
#     "tetra10"=> 10,
#     "hexahedron27"=> 27,
#     "wedge15"=> 15,
#     "wedge18"=> 18,
#     "pyramid13"=> 13,
#     "pyramid14"=> 14,
#     #
#     "line4"=> 4,
#     "triangle10"=> 10,
#     "quad16"=> 16,
#     "tetra20"=> 20,
#     "wedge40"=> 40,
#     "hexahedron64"=> 64,
#     #
#     "line5"=> 5,
#     "triangle15"=> 15,
#     "quad25"=> 25,
#     "tetra35"=> 35,
#     "wedge75"=> 75,
#     "hexahedron125"=> 125,
#     #
#     "line6"=> 6,
#     "triangle21"=> 21,
#     "quad36"=> 36,
#     "tetra56"=> 56,
#     "wedge126"=> 126,
#     "hexahedron216"=> 216,
#     #
#     "line7"=> 7,
#     "triangle28"=> 28,
#     "quad49"=> 49,
#     "tetra84"=> 84,
#     "wedge196"=> 196,
#     "hexahedron343"=> 343,
#     #
#     "line8"=> 8,
#     "triangle36"=> 36,
#     "quad64"=> 64,
#     "tetra120"=> 120,
#     "wedge288"=> 288,
#     "hexahedron512"=> 512,
#     #
#     "line9"=> 9,
#     "triangle45"=> 45,
#     "quad81"=> 81,
#     "tetra165"=> 165,
#     "wedge405"=> 405,
#     "hexahedron729"=> 729,
#     #
#     "line10"=> 10,
#     "triangle55"=> 55,
#     "quad100"=> 100,
#     "tetra220"=> 220,
#     "wedge550"=> 550,
#     "hexahedron1000"=> 1000,
#     "hexahedron1331"=> 1331,
#     #
#     "line11"=> 11,
#     "triangle66"=> 66,
#     "quad121"=> 121,
#     "tetra286"=> 286,
# )

# """
# I/O for Nastran bulk data.
# """
# function nasReader(filename::T) where {T<:AbstractString}

# nastran_to_meshio_type = Dict(  "CELAS1"=>"vertex",
#                                 "CBEAM"=>"line",
#                                 "CBUSH"=>"line",
#                                 "CBUSH1D"=>"line",
#                                 "CROD"=>"line",
#                                 "CGAP"=>"line",
#                                 "CBAR"=>"line",
#                                 "CTRIAR"=>"triangle",
#                                 "CTRIA3"=>"triangle",
#                                 "CTRAX6"=>"triangle6",
#                                 "CTRIAX6"=>"triangle6",
#                                 "CTRIA6"=>"triangle6",
#                                 "CQUADR"=>"quad",
#                                 "CSHEAR"=>"quad",
#                                 "CQUAD4"=>"quad",
#                                 "CQUAD8"=>"quad8",
#                                 "CQUAD9"=>"quad9",
#                                 "CTETRA"=>"tetra",
#                                 "CTETRA_"=>"tetra10",  # fictive
#                                 "CPYRAM"=>"pyramid",
#                                 "CPYRA"=>"pyramid",
#                                 "CPYRA_"=>"pyramid13",  # fictive
#                                 "CPENTA"=>"wedge",
#                                 "CPENTA_"=>"wedge15",  # fictive
#                                 "CHEXA"=>"hexahedron",
#                                 "CHEXA_"=>"hexahedron20",  # fictive
#                             )
# nastran_solid_types = ["CTETRA", "CPYRA", "CPENTA", "CHEXA"]
# meshio_to_nastran_type = Dict(v => k for (k, v) in nastran_to_meshio_type)


# function read(filename)
#     out = open(filename, "r") do f
#         read_buffer(f)
#     end
#     return out
# end


# function read_buffer(f)
#     # Skip until GRID
#     line =  readline(f)
#     while true
#         # nothing found
#         eof(f) && throw("GRID statement not found")
#         startswith(strip(line), "GRID") && break
#         line =  readline(f)
#     end

#     # Reading data
    
#      = []
#     points_id = []
#     cells = []
#     cells_id = []
#     cell = nothing
#     point_refs = []
#     cell_refs = []
#     cell_ref = nothing

#     function add_cell(keyword, cell, cell_ref)
#         cell_type = nastran_to_meshio_type[keyword]
#         cell = map(x -> parse(Int, x), cell)

#         # Treat 2nd order CTETRA, CPYRA, CPENTA, CHEXA elements
#         nastran_type = keyword
#         if length(cell) > num_nodes_per_cell[cell_type]
#             @assert cell_type in ["tetra", "pyramid", "wedge", "hexahedron"]
#             if cell_type == "tetra"
#                 cell_type = "tetra10"
#                 nastran_type = "CTETRA_"
#             elseif cell_type == "pyramid"
#                 cell_type = "pyramid13"
#                 nastran_type = "CPYRA_"
#             elseif cell_type == "wedge"
#                 cell_type = "wedge15"
#                 nastran_type = "CPENTA_"
#             elseif cell_type == "hexahedron"
#                 cell_type = "hexahedron20"
#                 nastran_type = "CHEXA_"
#             end
#         end

#         cell = _convert_to_vtk_ordering(cell, nastran_type)

#         if length(cells) > 0 && (cells[end][1] == cell_type)
#             push!(cells[end][2], cell)
#             push!(cells_id[end], cell_id)
#             if !isempty(cell_ref)
#                 push!(cell_refs, cell_ref)
#             end
#         else
#             push!(cells, (cell_type, [cell]))
#             push!(cells_id, [cell_id])
#             if !isnothing(cell_ref)
#                 push!(cell_refs, [cell_ref])
#             end
#         end
#     end # function

#     while true
#         # Blank lines or comments
#         if (length(line) < 4) || startswith(line, r"\$|//|#")
#             line = readline(f)
#             continue
#         else
#             break
#         end
#     end

#     while true
#         # End loop when ENDDATA detected
#         startswith(line, "ENDDATA") && break

#         # read line and merge with all continuation lines (starting with `+`)
#         chunks = _chunk_line(line)
#         while true
#             isempty(line) && throw(ReadError("Premature EOF"))
#             line = rstrip(line)
#             # Blank lines or comments
#             if length(line) < 4 || startswith(line, r"\$|//|#")
#                 continue
#             elseif line[1] == '+' || line[1] == '*'
#                 # skip the continuation chunk
#                 append!(chunks, _chunk_line(line))
#             else
#                 break
#             end
#             line = readline(f)
#         end

#         chunks  = [strip(chunk) for chunk in chunks]

#         keyword = chunks[1]

#         # Points
#         if keyword == "GRID"
#             point_id = parse(Int, chunks[2])
#             pref = strip(chunks[3])
#             length(pref) > 0 && append!(point_refs, parse(Int, pref))
#             append!(points_id, point_id)
#             append!(points, [_nastran_string_to_float(i) for i in chunks[4:6]])
#         elseif keyword == "GRID*"  # large field format: 8 + 16*4 + 8
#             point_id = parse(Int, chunks[2] * chunks[3])
#             pref = strip(chunks[4] * chunks[5])
#             if length(pref) > 0
#                 append!(point_refs, parse(Int, pref))
#                 append!(points_id, point_id)
#             end
#             chunks2 = _chunk_line(line)
#             line    =   readline(f)
#             append!(points, [_nastran_string_to_float(i * j) for (i, j) in [chunks[6:7], chunks[8:9], chunks2[2:3]]])
#         # CellBlock
#         elseif keyword in keys(nastran_to_meshio_type)
#             cell_id     =   parse(Int, chunks[2])
#             cell_refstring    =   strip(chunks[3])
#             cell_ref    =   (length(cell_refstring) > 0) ? parse(Int, cell_refstring) : nothing

#             cell = if keyword in ["CBAR", "CBEAM", "CBUSH", "CBUSH1D", "CGAP"]
#                 # Most Nastran 1D elements contain a third node (in the form of a node id or coordinates) to specify the local coordinate system:
#                 # https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/tdocExt/pdf/QRG.pdf
#                 # For example, a CBAR line can be
#                 # ```
#                 # CBAR          37               3       11.0     0.0     0.0
#                 # ```
#                 # where the last three floats specify the orientation vector.
#                 # This information is removed.
#                 chunks[4:5]
#             else
#                 chunks[4:end]
#             end

#             # remove empty chunks
#             cell = filter(!isempty, cell)

#             !isempty(cell) && add_cell(keyword, cell, cell_ref)
#         end # if
#     end # while

#     # Convert to arrays
#     points = Array(points)
#     points_id = Array{Int}(points_id)
#     for k, (c, cid) in enumerate(zip(cells, cells_id))
#         cells[k]    =   CellBlock(c[1], Array{Int}(c[2]))
#         cells_id[k] =   Array{Int}(cid)
#     end

#     # Convert to natural point ordering
#     # https://stackoverflow.com/questions/16992713/translate-every-element-in-numpy-array-according-to-key
#     points_id_dict = dict(zip(points_id, np.arange(len(points), dtype=int)))
#     points_id_get = np.vectorize(points_id_dict.__getitem__)
#     for k, c in enumerate(cells)
#         cells[k] = CellBlock(c.type, points_id_get(c.data))
#     end

#     # Construct the mesh object
#     mesh = Mesh(points, cells)
#     mesh.points_id = points_id
#     mesh.cells_id = cells_id
#     if len(point_refs) > 0:
#         mesh.point_data["nastran:ref"] = np.array(point_refs)
#     if len(cell_refs) > 0:
#         mesh.cell_data["nastran:ref"] = [np.array(i) for i in cell_refs]
#     return mesh


# # There are two basic categories of input data formats in NX Nastran:
# #
# #     "Free" format data, in which the data fields are simply separated by commas. This type of data is known as free field data.
# #     "Fixed" format data, in which your data must be aligned in columns of specific width. There are two subcategories of fixed format data that differ based on the size of the fixed column width:
# #         Small field format, in which a single line of data is divided into 10 fields that can contain eight characters each.
# #         Large field format, in which a single line of input is expanded into two lines The first and last fields on each line are eight columns wide, while the intermediate fields are sixteen columns wide. The large field format is useful when you need greater numerical accuracy.
# #
# # See: https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/tdocExt/pdf/User.pdf


# def write(filename, mesh, point_format="fixed-large", cell_format="fixed-small"):
#     if point_format == "free":
#         grid_fmt = "GRID,{:d},{:s},{:s},{:s},{:s}\n"
#         float_fmt = _float_to_nastran_string
#     elif point_format == "fixed-small":
#         grid_fmt = "GRID    {:<8d}{:<8s}{:>8s}{:>8s}{:>8s}\n"
#         float_fmt = _float_rstrip
#     elif point_format == "fixed-large":
#         grid_fmt = "GRID*   {:<16d}{:<16s}{:>16s}{:>16s}\n*       {:>16s}\n"
#         float_fmt = _float_to_nastran_string
#     else:
#         raise RuntimeError(f'unknown "{format}" format')

#     if cell_format == "free":
#         int_fmt, cell_info_fmt = "{:d}", "{:s},{:d},{:s},"
#         sjoin = ","
#     elif cell_format == "fixed-small":
#         int_fmt, cell_info_fmt = "{:<8d}", "{:<8s}{:<8d}{:<8s}"
#         sjoin, cchar = "", "+"
#         nipl1, nipl2 = 6, 14
#     elif cell_format == "fixed-large":
#         int_fmt, cell_info_fmt = "{:<16d}", "{:<8s}{:<16d}{:<16s}"
#         sjoin, cchar = "", "*"
#         nipl1, nipl2 = 2, 6
#     else:
#         raise RuntimeError(f'unknown "{format}" format')

#     if mesh.points.shape[1] == 2:
#         warn(
#             "Nastran requires 3D points, but 2D points given. "
#             "Appending 0 third component."
#         )
#         points = np.column_stack([mesh.points, np.zeros_like(mesh.points[:, 0])])
#     else:
#         points = mesh.points

#     with open_file(filename, "w") as f:
#         f.write(f"$ Nastran file written by meshio v{__version__}\n")
#         f.write("BEGIN BULK\n")

#         # Points
#         point_refs = mesh.point_data.get("nastran:ref", None)
#         for point_id, x in enumerate(points):
#             fx = [float_fmt(k) for k in x]
#             pref = str(point_refs[point_id]) if point_refs is not None else ""
#             f.write(grid_fmt.format(point_id + 1, pref, fx[0], fx[1], fx[2]))

#         # CellBlock
#         cell_id = 0
#         cell_refs = mesh.cell_data.get("nastran:ref", None)
#         for ict, cell_block in enumerate(mesh.cells):
#             cell_type = cell_block.type
#             cells = cell_block.data
#             nastran_type = meshio_to_nastran_type[cell_type]
#             if cell_format.endswith("-large"):
#                 nastran_type += "*"
#             if cell_refs is not None:
#                 cell_refs_t = cell_refs[ict]
#             else:
#                 cell_ref = ""
#             for ic, cell in enumerate(cells):
#                 if cell_refs is not None:
#                     cell_ref = str(int(cell_refs_t[ic]))
#                 cell_id += 1
#                 cell_info = cell_info_fmt.format(nastran_type, cell_id, cell_ref)
#                 cell1 = cell + 1
#                 cell1 = _convert_to_nastran_ordering(cell1, nastran_type)
#                 conn = sjoin.join(int_fmt.format(nid) for nid in cell1[:nipl1])
#                 if len(cell1) > nipl1:
#                     if cell_format == "free":
#                         cflag1 = cflag3 = ""
#                         cflag2 = cflag4 = "+,"
#                     else:
#                         cflag1 = cflag2 = f"{cchar}1{cell_id:<6x}"
#                         cflag3 = cflag4 = f"{cchar}2{cell_id:<6x}"
#                     f.write(cell_info + conn + cflag1 + "\n")
#                     conn = sjoin.join(int_fmt.format(nid) for nid in cell1[nipl1:nipl2])
#                     if len(cell1) > nipl2:
#                         f.write(cflag2 + conn + cflag3 + "\n")
#                         conn = sjoin.join(int_fmt.format(nid) for nid in cell1[nipl2:])
#                         f.write(cflag4 + conn + "\n")
#                     else:
#                         f.write(cflag2 + conn + "\n")
#                 else:
#                     f.write(cell_info + conn + "\n")

#         f.write("ENDDATA\n")


# def _float_rstrip(x, n=8):
#     return f"{x:f}".rstrip("0")[:n]


# def _float_to_nastran_string(value, length=16):
#     """
#     Return a value in NASTRAN scientific notation.
#     Examples:
#         1234.56789 --> "1.23456789+3"
#         -0.1234 --> "-1.234-1"
#         3.1415926535897932 --> "3.14159265359+0"
#     """
#     aux = length - 2
#     # sfmt = "{" + f":{length}s" + "}"
#     sfmt = "{" + ":s" + "}"
#     pv_fmt = "{" + f":{length}.{aux}e" + "}"

#     if value == 0.0:
#         return sfmt.format("0.")

#     python_value = pv_fmt.format(value)  # -1.e-2
#     svalue, sexponent = python_value.strip().split("e")
#     exponent = int(sexponent)  # removes 0s

#     sign = "-" if abs(value) < 1.0 else "+"

#     # the exponent will be added later...
#     sexp2 = str(exponent).strip("-+")
#     value2 = float(svalue)

#     # the plus 1 is for the sign
#     len_sexp = len(sexp2) + 1
#     leftover = length - len_sexp
#     leftover = leftover - 3 if value < 0 else leftover - 2
#     fmt = "{" + f":1.{leftover:d}f" + "}"

#     svalue3 = fmt.format(value2)
#     svalue4 = svalue3.strip("0")
#     field = sfmt.format(svalue4 + sign + sexp2)
#     return field


# function _nastran_string_to_float(string)
#     try
#         return parse(Float64, string)
#     catch ValueError
#         string = strip(string)
#         return float(string[1] * replace(string[1:end], "+"=>"e+", "-"=>"e-"))
#     end
# end

# function _chunk_line(line)
#     chunks = begin
#         if occursin(",", line)  # free format
#             split(line, ",")
#         else  # fixed format
#             CHUNK_SIZE = 8
#             @inbounds [line[i : CHUNK_SIZE + i - 1] for i in 1:CHUNK_SIZE:(length(line) - CHUNK_SIZE + 1)]
#         end
#     end # begin
#     # everything after the 9th chunk is ignored
#     return chunks[1:end]
# end


# function _convert_to_vtk_ordering(cell, nastran_type)
#     if nastran_type in ["CTRAX6", "CTRIAX6"]
#         cell = [cell[i] for i in [0, 2, 4, 1, 3, 5]]
#     elseif nastran_type == "CHEXA_"
#         cell = [ cell[i] for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15 ]]
#     elseif nastran_type == "CPENTA_"
#         cell = [cell[i] for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11]]
#     end
#     return cell
# end


# function _convert_to_nastran_ordering(cell, nastran_type)
#     if nastran_type in ["CTRAX6", "CTRIAX6"]
#         cell = [cell[i] for i in [0, 3, 1, 4, 2, 5]]
#     elseif nastran_type == "CHEXA_"
#         cell = [cell[i] for i in [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15]]
#     elseif nastran_type == "CPENTA_"
#         cell = [cell[i] for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11]]
#     end
#     return cell
# end


# register_format("nastran", [".bdf", ".fem", ".nas"], read, {"nastran": write})