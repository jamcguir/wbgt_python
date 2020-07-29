; This NCL script crops the elevation data to desired region, then regrids it to 2.5km to match RTMA & NDFD
;Created by Darrian Bertrand, 2019


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;Step 1: crop data
;Refer to regridding google doc if need help understanding process
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

f = addfile("/data/sco/wbgt_python/resources/elevation.nc4", "r")
elev = f->Band1
lat = f->lat
lon = f->lon

printVarSummary(elev)

elev_extract = elev({25.5:40.3277},{-88.5:-74.1883})   ;SE US
;elev_extract = elev({33.8189:39.5453},{-84.401:-75.2297})   ;all of NC & VA
;elev_extract = elev({33.8189:38.0664},{-84.401:-75.4563})  ;up to Charlottesville, VA

printVarSummary(elev_extract)

;---Write to netcdf

;---Open new file to write to 
;setfileoption("nc","preFill",False)
setfileoption("nc","defineMode",True)

;---Open a new NetCDF file to write to
  fout_name = "elevation_sercc.nc"
  ;fout_name = "elevation_crop_NC.nc"
  system("rm -f " + fout_name)
  fout = addfile(fout_name,"c")

;--Define dimensions
  nlat = 687    ;all of NC/VA
  nlon =  1100 
  ;nlat = 510   ;all of NC, to Charlottesville, VA
  ;nlon = 1072

;---Define the dimension names and their sizes on the file
  dims          = dimsizes(elev_extract)     ; nlat x nlev
  dim_names     = getvardims(elev_extract)   ; "lat", "lon"
  dim_unlimited = (/False,False/)
  filedimdef(fout,dim_names,dims,dim_unlimited)

;---Define each variable, its type, its dimension names, and its size
  filevardef(fout,"lat",typeof(lat),"lat")
  filevardef(fout,"lon",typeof(lon),"lon")
  filevardef(fout,"elevation",typeof(elev_extract),dim_names)

;---Define each variable's attributes.
  filevarattdef(fout,"lat",elev_extract&lat)
  filevarattdef(fout,"lon",elev_extract&lon)
  filevarattdef(fout,"elevation",elev_extract)

;---NOW write the variables to the file.
  fout->lat   = (/elev_extract&lat/)
  fout->lon   = (/elev_extract&lon/)
  fout->elev = elev_extract
  
  printVarSummary(elev_extract)
  ;exit()
  
;---Close file. Not necessary, but a good idea.
  delete(fout)


;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
; Step 2: Regrid using NCL's ESMF 
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

;--Data file containing source grid
    ;src_file = "/home/users/dmbertra/WBGT/elevation_crop.nc"
    src_file = "/data/sco/wbgt_python/elevation_sercc.nc"    ;file we created in Step 1
    sfile    = addfile(src_file,"r")
	src_lat = sfile->lat
	src_lon = sfile->lon
	
;---Get variable to regrid
    varname = "elev"                     
    var     = sfile->$varname$         ; Assumption is that "var"
                                       ; contains coordinate arrays.
    printVarSummary(var)               ; Make sure it has coord arrays.

;
; v2 - Read original data set 
;
    dst_file = "/data/sco/wbgt_python/input/rtma.nc"
    dfile    = addfile(dst_file,"r")
    dst_lat  = dfile->latitude                     ;;---Change (likely)
    dst_lon  = dfile->longitude                     ;;---Change (likely)


;---Set up regridding options
    Opt                   = True

;---"bilinear" is the default. "patch" and "conserve" are other options.
    Opt@InterpMethod      = "bilinear"        ;;---Change (maybe)

    Opt@WgtFileName       = "rect_to_curv.nc"

;
; These next two lines only needed if "var" doesn't
; contain coordinate arrays.
;
;   Opt@SrcGridLat        = sfile->lat        ;;--Change (maybe)
;   Opt@SrcGridLon        = sfile->lon        ;;--Change (maybe)

    Opt@SrcRegional       = True              ;;--Change (maybe)
    Opt@SrcInputFileName  = src_file          ; optional, but good idea
    Opt@SrcMask2D         = where(.not.ismissing(var),1,0) ; Necessary if has
                                                           ; missing values.

    Opt@DstGridLat        = dst_lat           ; destination grid
    Opt@DstGridLon        = dst_lon
    Opt@DstRegional       = True              ;;--Change (maybe)
    Opt@DstMask2D = where(.not.ismissing(dst_lat).and.\
                          .not.ismissing(dst_lon),1,0) ; Necessary if lat/lon
                                                      ; has missing values.

    Opt@ForceOverwrite    = True
    Opt@PrintTimings      = True
    Opt@Debug             = True

    var_regrid = ESMF_regrid(var,Opt)     ; Do the regridding

    printVarSummary(var_regrid)           ; Check that everything
    printMinMax(var_regrid,0)             ; looks okay.
	;exit()
;----------------------------------------------------------------------------
;-- Write to netcdf

;---Open new file to write to 
setfileoption("nc","defineMode",True)

;---Open a new NetCDF file to write to
  rgrdFileName = "elevation_regrid_sercc.nc"
  system("rm -f " + rgrdFileName)
  rgrd_nc = addfile(rgrdFileName,"c")

;---Create variable to hold global file attributes
      global = True
      copy_VarAtts(sfile, global)

      if (isatt(sfile,"title")) then
        global@TITLE = "REMAPPED: " + sfile@title
      end if
      global@remap         = "NCL: ESMF_regrid_with_weights (NCL version '" + \
                             get_ncl_version() + "')"
      global@remap_method  = "bilinear"
      global@creation_date = systemfunc("date")
    
      fileattdef( rgrd_nc, global )        ; copy global file attributes
      filedimdef(rgrd_nc,"TIME",-1,True)   ; force an unlimited dimension
;
; Write variables to file. Coordinate arrays will be written 
; automatically
;
      rgrd_nc->var = var_regrid
      rgrd_nc->lat = dst_lat
      rgrd_nc->lon = dst_lon
     
    
;---Get dimension sizes of regridded variable
      dims_regrid = dimsizes(var_regrid)


;end
