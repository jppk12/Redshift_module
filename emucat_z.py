#--------------------------------------------------#
# Author: Jahang Prathap P K
# Date: 31 January 2024
# Project: Redshift module for the EMUCAT pipeline
#--------------------------------------------------#
# Importing the necessary packages
#-----------------------------------#
import numpy as np, argparse
from astropy.table import Table, unique, vstack, join
import astropy.units as u
from astropy.coordinates import SkyCoord
from pyvo.dal import TAPService
import time
import pandas as pd
pd.pandas.set_option('display.max_columns', None)
#-----------------------------------#
# Defining an argument parser for 
# command line utility
#-----------------------------------#
def parse_args():
	'''
	Argument parser to take in arguments from 
	the command line
	'''
	parser=argparse.ArgumentParser(description="""To query tables in the
									EMU database and external ones to
									get the redshifts of the radio sources""")
	parser.add_argument("targets", help="""catalog of positions to query in 
                        votable format""")
	parser.add_argument("--upra", action='store', type=str, default='RAJ2000',
						help="The WISE RA column in the target list")
	parser.add_argument("--updec", action='store', type=str, default='DEJ2000',
						help="The WISE Dec column in the target list")
	parser.add_argument("--radio_ra", action='store', type=str, default='RA',
						help="The radio RA column in the target list")
	parser.add_argument("--radio_dec", action='store', type=str, default='DEC',
						help="The radio DEC column in the target list")
	parser.add_argument("--namecol",action='store',type=str,default='AllWISE',
						help="""ID column in targets. The primary key to look for 
						in EMUCAT""")
	parser.add_argument("--timeout",action='store',type=int,default=600,
						help="The query timeout limit")
	parser.add_argument("--searchrad",action='store',type=str,default='5arcsec',
						help="The crossmatch radius in the form xarcsec")
	parser.add_argument("--outdir",action='store',type=str,default='.',
						help="""The output directory to save the corssmatch results.
						Default is the current directory""")
						
	args=parser.parse_args()
	args.searchrad = u.Quantity(args.searchrad)
	args.outdir = args.outdir.removesuffix('/')
	
	return args
    
#----------------------------------------#
# Query the EMUCAT for Ned Taylor's
# redshift catalogue
#----------------------------------------#
# The functions here are written under the assumption that
# the input cataloge from EMUCat contains the WISE and radio
# coordinates. If there are some sources with no WISE counter-
# parts, there radio positions should still be available in the
# input catalogue. The search for redshifts are done in the order
# preference WISE coords first, if they are not available, then
# radio coords will be used.
#----------------------------------------#
def specz_ned(data,upra='RAJ2000',updec='DEJ2000',
            radio_ra='ra_deg_cont',radio_dec='dec_deg_cont',
            namecol='AllWISE',query_table='emucat.hemispec',wise_flag='wise_flag',
            searchrad=5*u.arcsec,qra='ra_deg',qdec='dec_deg',
            upload_name='upload_data', sep_col='angDist',finalcols=None):
    
    '''
    This function queries the EMUCat tap service to get the spectroscopic
    redhsifts for EMU sources. 
    '''
    t0=time.time()    
    data=Table.read(data)
    search_rad=str(searchrad.value)
    
    data.rename_column(name=upra, new_name=f'{upra}_wise')
    data.rename_column(name=updec, new_name=f'{updec}_wise')
    data.rename_column(name=radio_ra, new_name=f'{radio_ra}_radio')
    data.rename_column(name=radio_dec, new_name=f'{radio_dec}_radio')
    
    upra = f'{upra}_wise'
    updec = f'{updec}_wise'
    radio_ra = f'{radio_ra}_radio'
    radio_dec = f'{radio_dec}_radio'
    
    tab1 = data[data['AllWISE']!='']
    tab0 = data[data['AllWISE']=='']

    # Those with AllWISE coords
    if len(tab1) > 0:
    	query=f"""SELECT * FROM tap_upload.{upload_name} AS tup 
    		JOIN {query_table} AS t1
            ON 1=CONTAINS(POINT('ICRS',tup.{upra},tup.{updec}),
                            CIRCLE('ICRS',t1.{qra},t1.{qdec},{search_rad}/3600.)) 
            """
    	tap_emu = TAPService('https://emucat.aussrc.org/tap')
    	job = tap_emu.run_async(query,uploads={upload_name:tab1})
    
    	outdata1= job.to_table()
    	up_pos1 = SkyCoord(ra=outdata1[upra],dec=outdata1[updec])
    	q_pos1 = SkyCoord(ra=outdata1[qra],dec=outdata1[qdec])
    	angDist1 = up_pos1.separation(q_pos1).to(searchrad.unit)
    
    	outdata1[sep_col]=angDist1
    	outdata1[sep_col].unit=searchrad.unit
    else:
    	outdata1 = tab1
    	
    #Those without AllWISE coords
    if len(tab0) > 0:
    	query=f"""SELECT * FROM tap_upload.{upload_name} AS tup 
    		 JOIN {query_table} AS t1
            ON 1=CONTAINS(POINT('ICRS',tup.{radio_ra},tup.{radio_dec}),
                            CIRCLE('ICRS',t1.{qra},t1.{qdec},{search_rad}/3600.)) 
            """ 
    	job = tap_emu.run_async(query,uploads={upload_name:tab0})
    
    	outdata0= job.to_table()
    	up_pos0 = SkyCoord(ra=outdata0[radio_ra],dec=outdata0[radio_dec])
    	q_pos0 = SkyCoord(ra=outdata0[qra],dec=outdata0[qdec])
    	angDist0 = up_pos0.separation(q_pos0).to(searchrad.unit)
    
    	outdata0[sep_col]=angDist0
    	outdata0[sep_col].unit=searchrad.unit    
    else:
    	outdata0 = tab0
    
    # Stacking the tables together
    
    outdata = vstack([outdata1,outdata0])
    finalcols=[sep_col,namecol,upra,updec,radio_ra,radio_dec,
                'index','z_helio','z_quality','usez']
    outdata=outdata[finalcols]
    
    outdata.sort(sep_col)
    
    if len(outdata) <= 1:
    	outdata.rename_column(name=sep_col, new_name=f'{sep_col}_lit')
    	dt = time.time() - t0
    	print(f'''It took {np.round(dt,2)} s to find the hosts of {len(data)} sources from tap_emu''')
    	return outdata
    else:
        outdata = unique(outdata, radio_ra)
        outdata.sort(sep_col)
        outdata.rename_column(name=sep_col, new_name=f'{sep_col}_lit')
        dt = time.time() - t0
        print(f'''It took {np.round(dt,2)} s to find the hosts of {len(data)} sources from tap_emu''')
        return outdata
        
#----------------------------------------#
# Query the EMUCAT for Kenneth Duncan's
# LS DR8 photo-z estimates
#----------------------------------------#
def photoz_lsdr8(data, upra='RAJ2000',updec='DEJ2000',
                radio_ra='ra_deg_cont',radio_dec='dec_deg_cont',
                namecol='AllWISE',query_table='emucat.lsdr8_photoz',
                searchrad=5*u.arcsec,qra='ra',qdec='dec',upload_name='upload_data',
                sep_col='angDist',finalcols=None, wise_flag='wise_flag'):
    '''
    When possible, it is desirable to go for an IR detected radio source.
    The EMUCAT TAPService already has the emucat.lsdr8_photoz_nearest_allwise
    table which contains the Photo_z from Ken Duncan's LS DR8 and corresponding
    WISE detections. From a given target catalogue, the code below finds the WISE
    match first, get the photoz_id and assign the corresponding photo_z.
    '''
    t0=time.time()
    data=Table.read(data)
    search_rad=str(searchrad.value)
    
    data.rename_column(name=upra, new_name=f'{upra}_wise')
    data.rename_column(name=updec, new_name=f'{updec}_wise')
    data.rename_column(name=radio_ra, new_name=f'{radio_ra}_radio')
    data.rename_column(name=radio_dec, new_name=f'{radio_dec}_radio')
    
    upra = f'{upra}_wise'
    updec = f'{updec}_wise'
    radio_ra = f'{radio_ra}_radio'
    radio_dec = f'{radio_dec}_radio'
            
    tab1 = data[data['AllWISE']!='']
    tab0 = data[data['AllWISE']=='']
    
    # Those with AllWISE coords
    if len(tab1) > 0:
    	query=f"""
            SELECT * FROM tap_upload.{upload_name} AS tup JOIN {query_table} AS t1
            ON 1=CONTAINS(POINT('ICRS',tup.{upra},tup.{updec}),
                            CIRCLE('ICRS',t1.{qra},t1.{qdec},{search_rad}/3600.))    
            """
    	tap_emu = TAPService('https://emucat.aussrc.org/tap')
    	job = tap_emu.run_async(query,uploads={upload_name:tab1})
    	outdata1= job.to_table()
    
    	outdata1[qra].unit='deg'
    	outdata1[qdec].unit='deg'
    	up_pos1 = SkyCoord(ra=outdata1[upra],dec=outdata1[updec])
    	q_pos1 = SkyCoord(ra=outdata1[qra],dec=outdata1[qdec])
    	angDist1 = up_pos1.separation(q_pos1).to(searchrad.unit)
    
    	outdata1[sep_col]=angDist1
    	outdata1[sep_col].unit=searchrad.unit
    else:
    	outdata1 = tab1

    # Those without AllWISE coords
    if len(tab0) > 0:
    	query=f"""
            SELECT * FROM tap_upload.{upload_name} AS tup JOIN {query_table} AS t1
            ON 1=CONTAINS(POINT('ICRS',tup.{radio_ra},tup.{radio_dec}),
                            CIRCLE('ICRS',t1.{qra},t1.{qdec},{search_rad}/3600.))    
            """
    	tap_emu = TAPService('https://emucat.aussrc.org/tap')
    
    	job = tap_emu.run_async(query,uploads={upload_name:tab0})
    	outdata0= job.to_table()
    
    	outdata0[qra].unit='deg'
    	outdata0[qdec].unit='deg'
    	up_pos0 = SkyCoord(ra=outdata0[radio_ra],dec=outdata0[radio_dec])
    	q_pos0 = SkyCoord(ra=outdata0[qra],dec=outdata0[qdec])
    	angDist0 = up_pos0.separation(q_pos0).to(searchrad.unit)
    
    	outdata0[sep_col]=angDist0
    	outdata0[sep_col].unit=searchrad.unit
    else:
    	outdata0 = tab0

    # Stacking the tables together
    outdata = vstack([outdata1,outdata0])
    
    finalcols=[sep_col,namecol,upra,updec,radio_ra,radio_dec,'id','zphot','zphot_err']
    outdata = outdata[finalcols]
    outdata.sort(sep_col)
    if len(outdata) <= 1:
    	outdata.rename_column(name=sep_col,new_name=f'{sep_col}_lsdr8')
    	dt = time.time() - t0
    	print(f'''It took {np.round(dt,2)} s to find the hosts of {len(data)} sources from tap_emu''')
    	return outdata
    else:
        outdata = unique(outdata, radio_ra)
        outdata.sort(sep_col)
        outdata.rename_column(name=sep_col,new_name=f'{sep_col}_lsdr8')
        dt = time.time() - t0
        print(f'''It took {np.round(dt,2)} s to find the hosts of {len(data)} sources from tap_emu''')
        return outdata
        
#----------------------------------------#
# Query Vizier for MilliQuas redshifts
#----------------------------------------#
def milliquas_z(data,upra='RAJ2000',updec='DEJ2000',
                radio_ra='ra_deg_cont',radio_dec='dec_deg_cont',
                namecol='AllWISE',query_table='VII/294/catalog',
                searchrad=5*u.arcsec,qra='RAJ2000',wise_flag='wise_flag',
                qdec='DEJ2000',upload_name='upload_data',
                sep_col='angDist',finalcols=None):
    
    '''
    Queries the Million Quasar catalogue to get redshifts
    '''
    t0=time.time()
    data=Table.read(data)
    search_rad=str(searchrad.value)
        
    data.rename_column(name=upra, new_name=f'{upra}_wise')
    data.rename_column(name=updec, new_name=f'{updec}_wise')
    data.rename_column(name=radio_ra, new_name=f'{radio_ra}_radio')
    data.rename_column(name=radio_dec, new_name=f'{radio_dec}_radio')
    
    upra = f'{upra}_wise'
    updec = f'{updec}_wise'
    radio_ra = f'{radio_ra}_radio'
    radio_dec = f'{radio_dec}_radio'
    
    tab1 = data[data['AllWISE']!='']
    tab0 = data[data['AllWISE']=='']
    
    # Those with AllWISE coords
    if len(tab1) > 0:
    	query=f"""
            SELECT DISTANCE(POINT('ICRS',tup.{upra},tup.{updec}),
                            POINT('ICRS',t1.{qra},t1.{qdec})) AS angDist,
            tup.{namecol},tup.{upra},tup.{updec},tup.{radio_ra},tup.{radio_dec},
            t1.recno,t1.NAME,t1.z 
            FROM tap_upload.{upload_name} AS tup JOIN \"{query_table}\" 
            AS t1 ON 1=CONTAINS(POINT('ICRS',tup.{upra},tup.{updec}),
                            CIRCLE('ICRS',t1.{qra},t1.{qdec},{search_rad}/3600.))
          """
    	tap_vizier=TAPService('http://TAPVizieR.u-strasbg.fr/TAPVizieR/tap/')
    	job = tap_vizier.run_async(query,uploads={upload_name:tab1})
    
    	outdata1=job.to_table()
    	outdata1.sort(sep_col)
    else:
    	outdata1 = tab1

    # Those without AllWISE coords
    if len(tab0) > 0:
    	query=f"""
            SELECT DISTANCE(POINT('ICRS',tup.{radio_ra},tup.{radio_dec}),
                            POINT('ICRS',t1.{qra},t1.{qdec})) AS angDist,
            tup.{namecol},tup.{upra},tup.{updec},tup.{radio_ra},tup.{radio_dec},
            t1.recno,t1.NAME,t1.z 
            FROM tap_upload.{upload_name} AS tup JOIN \"{query_table}\" 
            AS t1 ON 1=CONTAINS(POINT('ICRS',tup.{radio_ra},tup.{radio_dec}),
                            CIRCLE('ICRS',t1.{qra},t1.{qdec},{search_rad}/3600.))
          """
    	tap_vizier=TAPService('http://TAPVizieR.u-strasbg.fr/TAPVizieR/tap/')
    	job = tap_vizier.run_async(query,uploads={upload_name:tab0})
    
    	outdata0=job.to_table()
    	outdata0.sort(sep_col)
    else:
    	outdata0 = tab0
	
    # Stacking the tables together
    outdata = vstack([outdata1,outdata0])
    outdata.sort(sep_col)
    
    if len(outdata) <= 1:
    	outdata.rename_column(name=sep_col, new_name=f'{sep_col}_mquas')
    	dt=time.time()-t0
    	print(f'''It took {np.round(dt,2)} s to find the hosts of {len(data)} sources from tap_vizier''')
    	return outdata
    else:
        outdata=unique(outdata,radio_ra)
        outdata.sort(sep_col)
        outdata.rename_column(name=sep_col, new_name=f'{sep_col}_mquas')
        dt=time.time()-t0
        print(f'''It took {np.round(dt,2)} s to find the hosts of {len(data)} sources from tap_vizier''')
        return outdata
        
#---------------------------------------#
# Implementing the command line utility
#---------------------------------------#
if __name__=='__main__':
    args = parse_args()
    data = args.targets
    
    specz_ned = specz_ned(data=data, upra=args.upra, updec=args.updec,
                          radio_ra=args.radio_ra, radio_dec=args.radio_dec,
                          namecol=args.namecol, searchrad=args.searchrad)
    photoz_lsdr8 = photoz_lsdr8(data=data, upra=args.upra, updec=args.updec,
                          radio_ra=args.radio_ra, radio_dec=args.radio_dec,
                          namecol=args.namecol, searchrad=args.searchrad)
    milliquas_z = milliquas_z(data=data, upra=args.upra, updec=args.updec,
                          radio_ra=args.radio_ra, radio_dec=args.radio_dec,
                          namecol=args.namecol, searchrad=args.searchrad)
    
    z_table = join(specz_ned['angDist_lit',f'{args.radio_ra}_radio',f'{args.radio_dec}_radio','index','z_helio','z_quality','usez'],milliquas_z['angDist_mquas',f'{args.radio_ra}_radio',f'{args.radio_dec}_radio','recno','Name','z'],join_type='outer')
    z_table = join(photoz_lsdr8, z_table, join_type = 'outer')
    z_table.write(f'{args.outdir}/z_table.xml',format='votable')
    
