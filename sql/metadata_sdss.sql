SELECT
  R.*,
  PO.ra as neighbor_ra, PO.dec as neighbor_dec,
  dbo.fDistanceArcMinEq(R.ra,R.dec,PO.ra,PO.dec) as sep_arcmin,
  PO.r as rmag, PO.err_r as err_rmag, 
  PO.type as sdsstype, 
  PZ.z AS neighbor_redshift, PZ.zErr AS neighbor_redshift_err, 
  PZ.absmagR, PZ.kcorrR,
  PZ.objid
  INTO mydb.bzcat_sdss
  FROM PhotoObj PO, mydb.all_blazars_neighbors R, PhotoZ PZ
WHERE R.matched_id = PO.objid
  AND PO.type = 3
  AND PZ.objid = R.matched_id
