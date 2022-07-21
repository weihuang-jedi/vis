;load "./create_netcdf_file_utils.ncl"
;
;start_write_time = get_cpu_time()

 wrkdir = "/scratch2/BMC/gsienkf/Wei.Huang/ufs/data/soca_letkf"
 typelist = (/"icec_ssmi_emc_2015120112", \
              "icefb_cryosat_esa_l2_2015120112", \
              "rads_adt_c2_2015120112", \
              "rads_adt_j2_2015120112", \
              "rads_adt_sa_2015120112", \
              "profile_wod_t_2015120112", \
              "profile_wod_s_2015120112", \
              "sst_avhrr19_l3u_nesdis_2015120112", \
              "sst_avhrrmta_l3u_nesdis_2015120112", \
              "sss_smap_jpl_L2_2015120112", \
              "sss_smos_esa_L2_2015120112"/)

 do nt=0,dimsizes(typelist)-1
    ftype = typelist(nt)

    print("ftype = " + ftype)

    cmdstr := "ls " + wrkdir + "/obsdir/output_" + ftype + "_*.nc"
    print("cmdstr=" + cmdstr)
    flst := systemfunc(cmdstr)
   ;print("flst=" + flst)
    fin := addfiles(flst, "r")

    vnames := getfilevarnames(fin[0])
   ;sqsort(vnames)

   ;ListSetType(fin, "cat")

   ;Open a new NetCDF file to write to
    setfileoption("nc","preFill",False)
    setfileoption("nc","defineMode",True)

   ;Open a new NetCDF file to write to
    fout_name = wrkdir + "/solver_input_" + ftype + ".nc"
    system("rm -f " + fout_name)
    fout = addfile(fout_name, "c")

    do nv=0,dimsizes(vnames)-1
       print("vname: " + vnames(nv))
       x := fin[:]->$vnames(nv)$
       fout->$vnames(nv)$ = x
    end do

   ;Close file. Not necessary, but a good idea.
    delete(fout)
 end do

;end_write_time = get_cpu_time()

