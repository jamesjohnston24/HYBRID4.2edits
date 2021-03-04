!======================================================================!
subroutine run_finish
use netcdf
use mpi
use shared
implicit none

character (len = 250) :: file_name ! Generic filename
integer :: k
integer :: ncid
integer :: varid
integer :: varid_kspp
integer :: varid_dbh,varid_hdbh
integer :: varid_cfoliage,varid_cfiner,varid_cstore
integer :: varid_nfoliage,varid_nbswood,varid_nheart,varid_navail
integer :: varid_nfiner
integer :: varid_lasa
integer :: varid_nitf
integer :: varid_hbc
integer :: varid_land_index
integer :: varid_cd,varid_dd,varid_bgs,varid_egs
integer :: varid_k_ind
integer :: varid_Cm,varid_Cu,varid_Cn,varid_Cv,varid_Ca,varid_Cs
integer :: varid_Nm,varid_Nu,varid_Nn,varid_Nv,varid_Na,varid_Ns
integer :: varid_Cpa
integer :: varid_Npa
integer :: varid_snmin
integer :: varid_nind
integer :: varid_snow,varid_soilw1,varid_soilw2,varid_soilw3
integer :: lon_dimid,lat_dimid
integer :: lon_varid,lat_varid
integer :: plot_dimid
integer :: plot_varid
integer :: ind_p_dimid,ind_t_dimid
integer :: ind_p_varid,ind_t_varid
integer :: land_dimid
integer :: land_varid
integer, dimension (1) :: dimids_one
integer, dimension (2) :: dimids_two
integer, dimension (3) :: dimids_three
integer, dimension (3) :: dimids_land

if (.NOT. (local)) then

!----------------------------------------------------------------------!
! Output the global LAI field.
!----------------------------------------------------------------------!
file_name = "LAI_grid.nc"
write (*, *) 'Writing to ', trim (file_name)
!----------------------------------------------------------------------!
! Create netCDF dataset and enter define mode.
!----------------------------------------------------------------------!
call check (nf90_create (trim (file_name), cmode = nf90_clobber, &
            ncid = ncid))
!----------------------------------------------------------------------!
! Define the dimensions.
!----------------------------------------------------------------------!
call check (nf90_def_dim (ncid, "longitude", nlon, lon_dimid))
call check (nf90_def_dim (ncid, "latitude" , nlat, lat_dimid))
!----------------------------------------------------------------------!
! Define coordinate variables.
!----------------------------------------------------------------------!
call check (nf90_def_var (ncid, "longitude", nf90_float, lon_dimid, &
            lon_varid))
call check (nf90_def_var (ncid, "latitude" , nf90_float, lat_dimid, &
            lat_varid))
dimids_two = (/ lon_dimid, lat_dimid /)
!----------------------------------------------------------------------!
! Assign units attributes to coordinate data.
!----------------------------------------------------------------------!
call check (nf90_put_att (ncid, lon_varid, "units", "degrees_east"))
call check (nf90_put_att (ncid, lat_varid, "units", "degrees_north"))
!----------------------------------------------------------------------!
! Define variable.
!----------------------------------------------------------------------!
call check (nf90_def_var (ncid, "LAI", nf90_float, &
            dimids_two, varid))
call check (nf90_put_att (ncid, varid, "units", "m2 m-2"))
call check (nf90_put_att (ncid, varid, "_FillValue", fillvalue))
!----------------------------------------------------------------------!
! End definitions.
!----------------------------------------------------------------------!
call check (nf90_enddef (ncid))
!----------------------------------------------------------------------!
! Write data.
!----------------------------------------------------------------------!
call check (nf90_put_var (ncid, lon_varid, lon))
call check (nf90_put_var (ncid, lat_varid, lat))
call check (nf90_put_var (ncid,     varid, LAI_grid))
!----------------------------------------------------------------------!
! Close file.
!----------------------------------------------------------------------!
call check (nf90_close (ncid))
!----------------------------------------------------------------------!
! Output global soil water fields.
!----------------------------------------------------------------------!
file_name = "Soil_water_grid.nc"
write (*, *) 'Writing to ', trim (file_name)
!----------------------------------------------------------------------!
! Create netCDF dataset and enter define mode.
!----------------------------------------------------------------------!
call check (nf90_create (trim (file_name), cmode = nf90_clobber, &
            ncid = ncid))
!----------------------------------------------------------------------!
! Define the dimensions.
!----------------------------------------------------------------------!
call check (nf90_def_dim (ncid, "longitude", nlon, lon_dimid))
call check (nf90_def_dim (ncid, "latitude" , nlat, lat_dimid))
!----------------------------------------------------------------------!
! Define coordinate variables.
!----------------------------------------------------------------------!
call check (nf90_def_var (ncid, "longitude", nf90_float, lon_dimid, &
            lon_varid))
call check (nf90_def_var (ncid, "latitude" , nf90_float, lat_dimid, &
            lat_varid))
dimids_two = (/ lon_dimid, lat_dimid /)
!----------------------------------------------------------------------!
! Assign units attributes to coordinate data.
!----------------------------------------------------------------------!
call check (nf90_put_att (ncid, lon_varid, "units", "degrees_east"))
call check (nf90_put_att (ncid, lat_varid, "units", "degrees_north"))
!----------------------------------------------------------------------!
! Define variables.
!----------------------------------------------------------------------!
call check (nf90_def_var (ncid, "Soil_water1", nf90_float, &
            dimids_two, varid1))
call check (nf90_put_att (ncid, varid1, "units", "m"))
call check (nf90_put_att (ncid, varid1, "_FillValue", fillvalue))
call check (nf90_def_var (ncid, "Soil_water2", nf90_float, &
            dimids_two, varid2))
call check (nf90_put_att (ncid, varid2, "units", "m"))
call check (nf90_put_att (ncid, varid2, "_FillValue", fillvalue))
call check (nf90_def_var (ncid, "Soil_water3", nf90_float, &
            dimids_two, varid3))
call check (nf90_put_att (ncid, varid3, "units", "m"))
call check (nf90_put_att (ncid, varid3, "_FillValue", fillvalue))
!----------------------------------------------------------------------!
! End definitions.
!----------------------------------------------------------------------!
call check (nf90_enddef (ncid))
!----------------------------------------------------------------------!
! Write data.
!----------------------------------------------------------------------!
call check (nf90_put_var (ncid, lon_varid, lon))
call check (nf90_put_var (ncid, lat_varid, lat))
call check (nf90_put_var (ncid,    varid1, soilw1))
call check (nf90_put_var (ncid,    varid2, soilw2))
call check (nf90_put_var (ncid,    varid3, soilw3))
!----------------------------------------------------------------------!
! Close file.
!----------------------------------------------------------------------!
call check (nf90_close (ncid))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Output the global NPP field.
!----------------------------------------------------------------------!
file_name = "NPP_grid.nc"
write (*, *) 'Writing to ', trim (file_name)
!----------------------------------------------------------------------!
! Create netCDF dataset and enter define mode.
!----------------------------------------------------------------------!
call check (nf90_create (trim (file_name), cmode = nf90_clobber, &
            ncid = ncid))
!----------------------------------------------------------------------!
! Define the dimensions.
!----------------------------------------------------------------------!
call check (nf90_def_dim (ncid, "longitude", nlon, lon_dimid))
call check (nf90_def_dim (ncid, "latitude" , nlat, lat_dimid))
!----------------------------------------------------------------------!
! Define coordinate variables.
!----------------------------------------------------------------------!
call check (nf90_def_var (ncid, "longitude", nf90_float, lon_dimid, &
            lon_varid))
call check (nf90_def_var (ncid, "latitude" , nf90_float, lat_dimid, &
            lat_varid))
dimids_two = (/ lon_dimid, lat_dimid /)
!----------------------------------------------------------------------!
! Assign units attributes to coordinate data.
!----------------------------------------------------------------------!
call check (nf90_put_att (ncid, lon_varid, "units", "degrees_east"))
call check (nf90_put_att (ncid, lat_varid, "units", "degrees_north"))
!----------------------------------------------------------------------!
! Define variable.
!----------------------------------------------------------------------!
call check (nf90_def_var (ncid, "NPP", nf90_float, &
            dimids_two, varid))
call check (nf90_put_att (ncid, varid, "units", "kg[C] m-2 yr-1"))
call check (nf90_put_att (ncid, varid, "_FillValue", fillvalue))
!----------------------------------------------------------------------!
! End definitions.
!----------------------------------------------------------------------!
call check (nf90_enddef (ncid))
!----------------------------------------------------------------------!
! Write data.
!----------------------------------------------------------------------!
call check (nf90_put_var (ncid, lon_varid, lon))
call check (nf90_put_var (ncid, lat_varid, lat))
call check (nf90_put_var (ncid,     varid, NPP_grid))
!----------------------------------------------------------------------!
! Close file.
!----------------------------------------------------------------------!
call check (nf90_close (ncid))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Output the global Cv field.
!----------------------------------------------------------------------!
file_name = "Cv_grid.nc"
write (*, *) 'Writing to ', trim (file_name)
!----------------------------------------------------------------------!
! Create netCDF dataset and enter define mode.
!----------------------------------------------------------------------!
call check (nf90_create (trim (file_name), cmode = nf90_clobber, &
            ncid = ncid))
!----------------------------------------------------------------------!
! Define the dimensions.
!----------------------------------------------------------------------!
call check (nf90_def_dim (ncid, "longitude", nlon, lon_dimid))
call check (nf90_def_dim (ncid, "latitude" , nlat, lat_dimid))
!----------------------------------------------------------------------!
! Define coordinate variables.
!----------------------------------------------------------------------!
call check (nf90_def_var (ncid, "longitude", nf90_float, lon_dimid, &
            lon_varid))
call check (nf90_def_var (ncid, "latitude" , nf90_float, lat_dimid, &
            lat_varid))
dimids_two = (/ lon_dimid, lat_dimid /)
!----------------------------------------------------------------------!
! Assign units attributes to coordinate data.
!----------------------------------------------------------------------!
call check (nf90_put_att (ncid, lon_varid, "units", "degrees_east"))
call check (nf90_put_att (ncid, lat_varid, "units", "degrees_north"))
!----------------------------------------------------------------------!
! Define variable.
!----------------------------------------------------------------------!
call check (nf90_def_var (ncid, "Cv", nf90_float, &
            dimids_two, varid))
call check (nf90_put_att (ncid, varid, "units", "kg[C] m-2"))
call check (nf90_put_att (ncid, varid, "_FillValue", fillvalue))
!----------------------------------------------------------------------!
! End definitions.
!----------------------------------------------------------------------!
call check (nf90_enddef (ncid))
!----------------------------------------------------------------------!
! Write data.
!----------------------------------------------------------------------!
call check (nf90_put_var (ncid, lon_varid, lon))
call check (nf90_put_var (ncid, lat_varid, lat))
call check (nf90_put_var (ncid,     varid, Cv_grid))
!----------------------------------------------------------------------!
! Close file.
!----------------------------------------------------------------------!
call check (nf90_close (ncid))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Output the global Cs field.
!----------------------------------------------------------------------!
file_name = "Cs_grid.nc"
write (*, *) 'Writing to ', trim (file_name)
!----------------------------------------------------------------------!
! Create netCDF dataset and enter define mode.
!----------------------------------------------------------------------!
call check (nf90_create (trim (file_name), cmode = nf90_clobber, &
            ncid = ncid))
!----------------------------------------------------------------------!
! Define the dimensions.
!----------------------------------------------------------------------!
call check (nf90_def_dim (ncid, "longitude", nlon, lon_dimid))
call check (nf90_def_dim (ncid, "latitude" , nlat, lat_dimid))
!----------------------------------------------------------------------!
! Define coordinate variables.
!----------------------------------------------------------------------!
call check (nf90_def_var (ncid, "longitude", nf90_float, lon_dimid, &
            lon_varid))
call check (nf90_def_var (ncid, "latitude" , nf90_float, lat_dimid, &
            lat_varid))
dimids_two = (/ lon_dimid, lat_dimid /)
!----------------------------------------------------------------------!
! Assign units attributes to coordinate data.
!----------------------------------------------------------------------!
call check (nf90_put_att (ncid, lon_varid, "units", "degrees_east"))
call check (nf90_put_att (ncid, lat_varid, "units", "degrees_north"))
!----------------------------------------------------------------------!
! Define variable.
!----------------------------------------------------------------------!
call check (nf90_def_var (ncid, "Cs", nf90_float, &
            dimids_two, varid))
call check (nf90_put_att (ncid, varid, "units", "kg[C] m-2"))
call check (nf90_put_att (ncid, varid, "_FillValue", fillvalue))
!----------------------------------------------------------------------!
! End definitions.
!----------------------------------------------------------------------!
call check (nf90_enddef (ncid))
!----------------------------------------------------------------------!
! Write data.
!----------------------------------------------------------------------!
call check (nf90_put_var (ncid, lon_varid, lon))
call check (nf90_put_var (ncid, lat_varid, lat))
call check (nf90_put_var (ncid,     varid, Cs_grid))
!----------------------------------------------------------------------!
! Close file.
!----------------------------------------------------------------------!
call check (nf90_close (ncid))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Output the global Cv_GPT field.
!----------------------------------------------------------------------!
file_name = "Cv_GPT_grid.nc"
write (*, *) 'Writing to ', trim (file_name)
!----------------------------------------------------------------------!
! Create netCDF dataset and enter define mode.
!----------------------------------------------------------------------!
call check (nf90_create (trim (file_name), cmode = nf90_clobber, &
            ncid = ncid))
!----------------------------------------------------------------------!
! Define the dimensions.
!----------------------------------------------------------------------!
call check (nf90_def_dim (ncid, "longitude", nlon, lon_dimid))
call check (nf90_def_dim (ncid, "latitude" , nlat, lat_dimid))
!----------------------------------------------------------------------!
! Define coordinate variables.
!----------------------------------------------------------------------!
call check (nf90_def_var (ncid, "longitude", nf90_float, lon_dimid, &
            lon_varid))
call check (nf90_def_var (ncid, "latitude" , nf90_float, lat_dimid, &
            lat_varid))
dimids_two = (/ lon_dimid, lat_dimid /)
!----------------------------------------------------------------------!
! Assign units attributes to coordinate data.
!----------------------------------------------------------------------!
call check (nf90_put_att (ncid, lon_varid, "units", "degrees_east"))
call check (nf90_put_att (ncid, lat_varid, "units", "degrees_north"))
!----------------------------------------------------------------------!
! Define variables.
!----------------------------------------------------------------------!
call check (nf90_def_var (ncid, "Cv_C3GR", nf90_float, &
            dimids_two, varid_C3GR))
call check (nf90_put_att (ncid, varid_C3GR, "units", "kg[C] m-2"))
call check (nf90_put_att (ncid, varid_C3GR, "_FillValue", fillvalue))

call check (nf90_def_var (ncid, "Cv_C4GR", nf90_float, &
            dimids_two, varid_C4GR))
call check (nf90_put_att (ncid, varid_C4GR, "units", "kg[C] m-2"))
call check (nf90_put_att (ncid, varid_C4GR, "_FillValue", fillvalue))

call check (nf90_def_var (ncid, "Cv_BRCD", nf90_float, &
            dimids_two, varid_BRCD))
call check (nf90_put_att (ncid, varid_BRCD, "units", "kg[C] m-2"))
call check (nf90_put_att (ncid, varid_BRCD, "_FillValue", fillvalue))

call check (nf90_def_var (ncid, "Cv_NLEV", nf90_float, &
            dimids_two, varid_NLEV))
call check (nf90_put_att (ncid, varid_NLEV, "units", "kg[C] m-2"))
call check (nf90_put_att (ncid, varid_NLEV, "_FillValue", fillvalue))
!----------------------------------------------------------------------!
! End definitions.
!----------------------------------------------------------------------!
call check (nf90_enddef (ncid))
!----------------------------------------------------------------------!
! Write data.
!----------------------------------------------------------------------!
call check (nf90_put_var (ncid, lon_varid, lon))
call check (nf90_put_var (ncid, lat_varid, lat))
call check (nf90_put_var (ncid,     varid_C3GR, Cv_C3GR_grid))
call check (nf90_put_var (ncid,     varid_C4GR, Cv_C4GR_grid))
call check (nf90_put_var (ncid,     varid_BRCD, Cv_BRCD_grid))
call check (nf90_put_var (ncid,     varid_NLEV, Cv_NLEV_grid))
!----------------------------------------------------------------------!
! Close file.
!----------------------------------------------------------------------!
call check (nf90_close (ncid))
!----------------------------------------------------------------------!

end if ! .NOT. (local)

!----------------------------------------------------------------------!
! Output state variables to restart file if required.
!----------------------------------------------------------------------!
if (rsf_out) then
 !---------------------------------------------------------------------!
 ! For coordinates in restart file.
 ! 4 byte integer can go up to 2,147,483,647 , I think.
 ! 10 plots over 67420 grid-boxes would allow 3185 individuals/plot.
 !---------------------------------------------------------------------!
 allocate (ind(nind_total))
 do k = 1, nind_total
  ind (k) = k
 end do
 !----------------------------------------------------------------------!

 !----------------------------------------------------------------------!
 file_name = &
            "/home/adf10/rds/rds-mb425-geogscratch/adf10/RSF/rsf_out.nc"
 write (*, *) 'Writing to ', trim (file_name)
 !---------------------------------------------------------------------!
 ! Create netCDF dataset and enter define mode.
 ! nf90_64bit_offset required because file is large. Found rather
 ! randomly. Allows it to work when added lasa to output.
 !---------------------------------------------------------------------!
 call check (nf90_create (trim (file_name), cmode = nf90_64bit_offset, &
             ncid = ncid))
 !---------------------------------------------------------------------!
 ! Define the dimensions.
 !---------------------------------------------------------------------!
 call check (nf90_def_dim (ncid, "longitude", nlon      , lon_dimid  ))
 call check (nf90_def_dim (ncid, "latitude" , nlat      , lat_dimid  ))
 call check (nf90_def_dim (ncid, "plot"     , nplots    , plot_dimid ))
 call check (nf90_def_dim (ncid, "ind_plot" , nind_max  , ind_p_dimid))
 call check (nf90_def_dim (ncid, "ind_total", nind_total, ind_t_dimid))
 call check (nf90_def_dim (ncid, "iland"    , nland_max , land_dimid ))
 !---------------------------------------------------------------------!
 ! Define coordinate variables.
 !---------------------------------------------------------------------!
 call check (nf90_def_var (ncid, "longitude", nf90_float, lon_dimid , &
             lon_varid ))
 call check (nf90_def_var (ncid, "latitude" , nf90_float, lat_dimid , &
             lat_varid ))
 call check (nf90_def_var (ncid, "plot"     , nf90_int  , plot_dimid, &
             plot_varid))
 call check (nf90_def_var (ncid, "ind_plot" , nf90_int  , ind_p_dimid, &
             ind_p_varid))
 call check (nf90_def_var (ncid, "ind_total", nf90_int  , ind_t_dimid, &
             ind_t_varid))
 call check (nf90_def_var (ncid, "iland"    , nf90_int  , land_dimid , &
             land_varid))
 !---------------------------------------------------------------------!
 dimids_one   = (/ ind_t_dimid /)
 dimids_two   = (/ lon_dimid, lat_dimid /)
 dimids_three = (/ lon_dimid, lat_dimid, plot_dimid /)
 dimids_land  = (/ land_dimid, plot_dimid, ind_p_dimid /)
 !---------------------------------------------------------------------!
 ! Assign units attributes to coordinate data.
 !--------------------------------------------------------------------!
 call check (nf90_put_att (ncid, lon_varid, "units", "degrees_east"))
 call check (nf90_put_att (ncid, lat_varid, "units", "degrees_north"))
 !---------------------------------------------------------------------!
 ! Define variables.
 !---------------------------------------------------------------------!
 call check (nf90_def_var (ncid, "Land_index", &
  nf90_int, dimids_two, varid_land_index))
 call check (nf90_def_var (ncid, "Chilling-days", &
  nf90_int, dimids_two, varid_cd))
 call check (nf90_def_var (ncid, "Degree-days", &
  nf90_float, dimids_two, varid_dd))
 call check (nf90_def_var (ncid, "Beginning_of_growing_season", &
  nf90_int, dimids_two, varid_bgs))
 call check (nf90_def_var (ncid, "End_of_growing_season", &
  nf90_int, dimids_two, varid_egs))
 call check (nf90_def_var (ncid, "Individual_index", &
  nf90_int, dimids_land, varid_k_ind))
 call check (nf90_def_var (ncid, "Surface_metabolic_organic_C", &
  nf90_float, dimids_three, varid_Cm))
 call check (nf90_def_var (ncid, "Surface_structural_organic_C", &
  nf90_float, dimids_three, varid_Cu))
 call check (nf90_def_var (ncid, "Soil_metabolic_organic_C", &
  nf90_float, dimids_three, varid_Cn))
 call check (nf90_def_var (ncid, "Soil_structural_organic_C", &
  nf90_float, dimids_three, varid_Cv))
 call check (nf90_def_var (ncid, "Active_organic_C", &
  nf90_float, dimids_three, varid_Ca))
 call check (nf90_def_var (ncid, "Slow_organic_C", &
  nf90_float, dimids_three, varid_Cs))
 call check (nf90_def_var (ncid, "Passive_soil_organic_C", &
  nf90_float, dimids_three, varid_Cpa))
 call check (nf90_def_var (ncid, "Surface_metabolic_organic_N", &
  nf90_float, dimids_three, varid_Nm))
 call check (nf90_def_var (ncid, "Surface_structural_organic_N", &
  nf90_float, dimids_three, varid_Nu))
 call check (nf90_def_var (ncid, "Soil_metabolic_organic_N", &
  nf90_float, dimids_three, varid_Nn))
 call check (nf90_def_var (ncid, "Soil_structural_organic_N", &
  nf90_float, dimids_three, varid_Nv))
 call check (nf90_def_var (ncid, "Active_organic_N", &
  nf90_float, dimids_three, varid_Na))
 call check (nf90_def_var (ncid, "Slow_organic_N", &
  nf90_float, dimids_three, varid_Ns))
 call check (nf90_def_var (ncid, "Passive_soil_organic_N", &
  nf90_float, dimids_three, varid_Npa))
 call check (nf90_def_var (ncid, "Mineral_N", &
  nf90_float, dimids_three, varid_snmin))
 call check (nf90_def_var (ncid, "Number_individuals_per_plot", &
  nf90_int, dimids_three, varid_nind))
 call check (nf90_def_var (ncid, "Snowpack", &
  nf90_float, dimids_three, varid_snow))
 call check (nf90_def_var (ncid, "Soil_water_layer_1", &
  nf90_float, dimids_three, varid_soilw1))
 call check (nf90_def_var (ncid, "Soil_water_layer_2", &
  nf90_float, dimids_three, varid_soilw2))
 call check (nf90_def_var (ncid, "Soil_water_layer_3", &
  nf90_float, dimids_three, varid_soilw3))
 call check (nf90_def_var (ncid, "Species_number", &
  nf90_int, dimids_one, varid_kspp))
 call check (nf90_def_var (ncid, "DBH", &
  nf90_float, dimids_one, varid_dbh))
 call check (nf90_def_var (ncid, "Heartwood_DBH", &
  nf90_float, dimids_one, varid_hdbh))
 call check (nf90_def_var (ncid, "Foliage_C", &
  nf90_float, dimids_one, varid_cfoliage))
 call check (nf90_def_var (ncid, "Fine_root_C", &
  nf90_float, dimids_one, varid_cfiner))
 call check (nf90_def_var (ncid, "Storage_C", &
  nf90_float, dimids_one, varid_cstore))
 call check (nf90_def_var (ncid, "Foliage_N", &
  nf90_float, dimids_one, varid_nfoliage))
 call check (nf90_def_var (ncid, "Non-heartwood_structure_N", &
  nf90_float, dimids_one, varid_nbswood))
 call check (nf90_def_var (ncid, "Heartwood_N", &
  nf90_float, dimids_one, varid_nheart))
 call check (nf90_def_var (ncid, "Storage_N", &
  nf90_float, dimids_one, varid_navail))
 call check (nf90_def_var (ncid, "Fine_root_N", &
  nf90_float, dimids_one, varid_nfiner))
 call check (nf90_def_var (ncid, "Foliage_sapwood_area_ratio", &
  nf90_float, dimids_one, varid_lasa))
 call check (nf90_def_var (ncid, "Foliage_cold_damage", &
  nf90_float, dimids_one, varid_nitf))
 call check (nf90_def_var (ncid, "Height_to_base_of_crown", &
  nf90_int  , dimids_one, varid_hbc))
 !---------------------------------------------------------------------!
 call check (nf90_put_att (ncid, varid_cd , "units", "days"))
 call check (nf90_put_att (ncid, varid_dd , "units", "oC"))
 call check (nf90_put_att (ncid, varid_bgs, "units", "DOY"))
 call check (nf90_put_att (ncid, varid_egs, "units", "DOY"))
 call check (nf90_put_att (ncid, varid_Cm , "units", "kg[C] m-2"))
 call check (nf90_put_att (ncid, varid_Cu , "units", "kg[C] m-2"))
 call check (nf90_put_att (ncid, varid_Cn , "units", "kg[C] m-2"))
 call check (nf90_put_att (ncid, varid_Cv , "units", "kg[C] m-2"))
 call check (nf90_put_att (ncid, varid_Ca , "units", "kg[C] m-2"))
 call check (nf90_put_att (ncid, varid_Cs , "units", "kg[C] m-2"))
 call check (nf90_put_att (ncid, varid_Cpa, "units", "kg[C] m-2"))
 call check (nf90_put_att (ncid, varid_Nm , "units", "kg[N] m-2"))
 call check (nf90_put_att (ncid, varid_Nu , "units", "kg[N] m-2"))
 call check (nf90_put_att (ncid, varid_Nn , "units", "kg[N] m-2"))
 call check (nf90_put_att (ncid, varid_Nv , "units", "kg[N] m-2"))
 call check (nf90_put_att (ncid, varid_Na , "units", "kg[N] m-2"))
 call check (nf90_put_att (ncid, varid_Ns , "units", "kg[N] m-2"))
 call check (nf90_put_att (ncid, varid_Npa, "units", "kg[N] m-2"))
 call check (nf90_put_att (ncid, varid_snmin, "units", "kg[N] m-2"))
 call check (nf90_put_att (ncid, varid_nind , "units", "n"))
 call check (nf90_put_att (ncid, varid_snow , "units", "m"))
 call check (nf90_put_att (ncid, varid_soilw1  , "units", "m"))
 call check (nf90_put_att (ncid, varid_soilw2  , "units", "m"))
 call check (nf90_put_att (ncid, varid_soilw3  , "units", "m"))
 call check (nf90_put_att (ncid, varid_kspp    , "units", "n"))
 call check (nf90_put_att (ncid, varid_dbh     , "units", "m"))
 call check (nf90_put_att (ncid, varid_hdbh    , "units", "m"))
 call check (nf90_put_att (ncid, varid_cfoliage, "units", "kg[C] ind-1"))
 call check (nf90_put_att (ncid, varid_cfiner  , "units", "kg[C] ind-1"))
 call check (nf90_put_att (ncid, varid_cstore  , "units", "kg[C] ind-1"))
 call check (nf90_put_att (ncid, varid_nfoliage, "units", "kg[N] ind-1"))
 call check (nf90_put_att (ncid, varid_nbswood , "units", "kg[N] ind-1"))
 call check (nf90_put_att (ncid, varid_nheart  , "units", "kg[N] ind-1"))
 call check (nf90_put_att (ncid, varid_navail  , "units", "kg[N] ind-1"))
 call check (nf90_put_att (ncid, varid_nfiner  , "units", "kg[N] ind-1"))
 call check (nf90_put_att (ncid, varid_lasa    , "units", "m2 m-2"))
 call check (nf90_put_att (ncid, varid_nitf    , "units", "1-fraction"))
 call check (nf90_put_att (ncid, varid_hbc     , "units", "m"))
 !---------------------------------------------------------------------!
 ! End definitions.
 !---------------------------------------------------------------------!
 call check (nf90_enddef (ncid))
 !---------------------------------------------------------------------!
 ! Write data.
 !---------------------------------------------------------------------!
 call check (nf90_put_var (ncid, lon_varid  ,   lon))
 call check (nf90_put_var (ncid, lat_varid  ,   lat))
 call check (nf90_put_var (ncid, plot_varid ,  plot))
 call check (nf90_put_var (ncid, ind_t_varid,   ind))
 call check (nf90_put_var (ncid, varid_land_index, land_index))
 call check (nf90_put_var (ncid, varid_cd, cd))
 call check (nf90_put_var (ncid, varid_dd, dd))
 call check (nf90_put_var (ncid, varid_bgs, bgs))
 call check (nf90_put_var (ncid, varid_egs, egs))
 call check (nf90_put_var (ncid, varid_k_ind, k_ind))
 call check (nf90_put_var (ncid, varid_Cm , Cm))
 call check (nf90_put_var (ncid, varid_Cu , Cu))
 call check (nf90_put_var (ncid, varid_Cn , Cn))
 call check (nf90_put_var (ncid, varid_Cv , Cv))
 call check (nf90_put_var (ncid, varid_Ca , Ca))
 call check (nf90_put_var (ncid, varid_Cs , Cs))
 call check (nf90_put_var (ncid, varid_Cpa, Cpa))
 call check (nf90_put_var (ncid, varid_Nm , Nm))
 call check (nf90_put_var (ncid, varid_Nu , Nu))
 call check (nf90_put_var (ncid, varid_Nn , Nn))
 call check (nf90_put_var (ncid, varid_Nv , Nv))
 call check (nf90_put_var (ncid, varid_Na , Na))
 call check (nf90_put_var (ncid, varid_Ns , Ns))
 call check (nf90_put_var (ncid, varid_Npa, Npa))
 call check (nf90_put_var (ncid, varid_snmin, snmin))
 call check (nf90_put_var (ncid, varid_nind, nind))
 call check (nf90_put_var (ncid, varid_snow, snow))
 call check (nf90_put_var (ncid, varid_soilw1, soilw1))
 call check (nf90_put_var (ncid, varid_soilw2, soilw2))
 call check (nf90_put_var (ncid, varid_soilw3, soilw3))
 call check (nf90_put_var (ncid, varid_kspp  , kspp  ))
 call check (nf90_put_var (ncid, varid_dbh     , dbh))
 call check (nf90_put_var (ncid, varid_hdbh    , hdbh))
 call check (nf90_put_var (ncid, varid_cfoliage, cfoliage))
 call check (nf90_put_var (ncid, varid_cfiner  , cfiner))
 call check (nf90_put_var (ncid, varid_cstore  , cstore))
 call check (nf90_put_var (ncid, varid_nfoliage, nfoliage))
 call check (nf90_put_var (ncid, varid_nbswood , nbswood))
 call check (nf90_put_var (ncid, varid_nheart  , nheart))
 call check (nf90_put_var (ncid, varid_navail  , navail))
 call check (nf90_put_var (ncid, varid_nfiner  , nfiner))
 call check (nf90_put_var (ncid, varid_lasa    , lasa))
 call check (nf90_put_var (ncid, varid_nitf    , nitf))
 call check (nf90_put_var (ncid, varid_hbc     , hbc))
 !---------------------------------------------------------------------!
 ! Close file.
 !---------------------------------------------------------------------!
 call check (nf90_close (ncid))
 !---------------------------------------------------------------------!
end if

if ((local) .and. (wrclim)) close (20) ! Local climate output.
if (local) close (21) ! soil carbon output.
if (local) close (22) ! soil water output.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Close diagnostic output files.
!----------------------------------------------------------------------!
close (20) ! file='annual_global.dat'
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
contains
 subroutine check (status)

 integer, intent (in) :: status
 if (status /= nf90_noerr) then
  print *, trim (nf90_strerror(status))
  stop  "Stopped"
 end if
 end subroutine check
!----------------------------------------------------------------------!

end subroutine run_finish
!======================================================================!
