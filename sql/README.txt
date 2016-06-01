K. Willett (UMN) - 1 Jun 2016

There are two SQL queries run on CasJobs (http://skyserver.sdss.org/casjobs/)
used for studying blazar environments. 

1. neighbors_sdss.sql

    The first is neighbors_sdss.sql. This used to be accessible through a button
    called "Neighbors" on the MyDB page; as of Jun 2016, it seems to have been
    deprecated. The query should still run if fully executed, though. 
    
    The intent of this query is to take a list of blazars that have already
    been identified and preloaded from BZCAT5 (all_blazars) and search for nearby
    neighbors. I don't know of an easy way to automatically limit this to particular
    types of objects (such as galaxies), so this query searches for everything
    and we'll cut later. 

    The maximum search radius in this example is set at 10 arcmin, which will
    effectively set the upper redshift limit of our sample (since we want to
    keep some fixed physical radius to have enough area around each galaxy for
    a robust measure of environment). If it's much larger than this, CasJobs
    returns too many objects to practically manage the data (this already 
    returns about 6 million rows). Thankfully, we don't need to make it much bigger
    since the radius corresponding to 50 Mpc at 10 arcmin gets close to the
    magnitude limit of SDSS. 

    This query took about 90 minutes to run on CasJobs.

2. metadata_sdss.sql

    The second SQL query takes every matched object from the first and does two
    things: limits it to only galaxies (so we throw out stars and quasars)
    and retrieves metadata that'll be necessary for calculating B_gB. The
    critical things we need are magnitude and redshift, but a number of
    other photometric parameters are included in case we want to do
    further analysis later.  

    This reduces the full sample from about 1 million rows to 6 million rows. 

    This query took about 5 minutes to run on CasJobs.

Once both queries are done, we want to download the data from all_blazars_sdss
locally so we can calculate the clustering amplitude. The export routines on
CasJobs have been buggy and broken for some time (despite repeated complaints)
- in particular, it won't let you download it as a FITS binary table. To keep
things like the data type fixed, I usually download it as an "XML - Virtual
Observatory VOTABLE". I then open it in TOPCAT and re-save as a FITS binary,
which is what bgb.py will expect. It's also renamed: the file
blazars_unique_all_10arcmin_sdss.fits is the current input to bgb.py, 
and what you should be able to recreate following the above steps. 
