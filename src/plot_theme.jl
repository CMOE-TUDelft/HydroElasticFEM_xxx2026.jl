"""
Shared CairoMakie theme for all verification-suite figures
(same palette and styling as HydroElasticFEM_Colomes2023.jl).
"""

const VLFS_THEME = Theme(
  fontsize = 12,
  palette = (color = [:black, :royalblue, :firebrick, :forestgreen, :darkorange],),
  Axis = (
    xgridvisible = false,
    ygridvisible = false,
  ),
  Lines = (
    linewidth = 1.8,
  ),
  Scatter = (
    markersize = 8,
  ),
)
