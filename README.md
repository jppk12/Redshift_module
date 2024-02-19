This module contains a set of functions that query the EMUCat and VizieR
TAP Services for redshifts. Three functions are available now to cross-match 
a given input table with the 'emucat.hemispec', 'emucat.lsdr8_photoz', and
'VII/294/catalog' (the milliquasar redshift catalogue from VizieR.)

The functions expect a votable with the equatorial coordinates of the source of 
interest. The cross-match is done using the AllWISE positions of the sources,
if available, else, the radio positions will be used. 

The functions can be used from the command-line. An example execution is as 
follows:
For a file with name query_data.xml with the WISE coords RAJ2000 and DEJ2000
and radio coords as ra_deg_cont and dec_deg_cont,

python emucat_z.py query_data.xml --upra=RAJ2000 --updec=DEJ2000
 --radio_ra=ra_deg_cont --radio_dec=dec_deg_cont --searchrad=3arcsec 

The code outputs a votable that contains the separation from all three tables,
the AllWISE ID, the WISE and radio coords, the redshift and related parameters.

An example input file (query_data.xml) and an example output file 
(z_table.xml) is also made available.
