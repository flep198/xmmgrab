import requests
import pandas as pd
from datetime import datetime
import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
from tqdm import tqdm
from pandas.io.html import read_html
from selenium import webdriver
import contextlib
import csv, smtplib, ssl
import time
from astropy.time import Time
from astroquery.ipac.ned import Ned

#Telescopes to consider
effelsberg=EarthLocation(lat=50.5*u.deg, lon=6.9*u.deg, height=319*u.m)
eff_elev_lim=10*u.deg
parkes=EarthLocation(lat=-33*u.deg, lon=148.25*u.deg, height=324*u.m)
parkes_elev_lim=30*u.deg
fast=EarthLocation(lat=25.65*u.deg, lon=106.86*u.deg)
fast_elev_lim=60*u.deg

telescopes=[{"name": "Effelsberg", "lat": 50.5*u.deg, "lon": 6.9*u.deg, "height": 319*u.m, "elev_lim": 10*u.deg},
        {"name": "Parkes", "lat": -33*u.deg, "lon": 148.25*u.deg, "height": 324*u.m, "elev_lim": 30*u.deg},
        {"name": "FAST", "lat": 25.65*u.deg, "lon": 106.86*u.deg, "elev_lim": 60*u.deg},
        {"name": "MeerKAT", "lat": -30.711*u.deg,"lon": 21.44*u.deg,"elev_lim": 15*u.deg}]


def checkObs():
    output=[]
    
    #list of target names to be observed
    target_file=open('xmm_targets.csv')
    target_names=[]
    for row in target_file:
        target_names=np.append(target_names,row)

    url = 'https://xmm-tools.cosmos.esa.int/external/xmm_sched/short_term_schedule.php'
    html = requests.get(url).content
    df_list = pd.read_html(html,header =0)
    df = pd.DataFrame(df_list[-1])

    for ind,name in enumerate(df["Target Name"]):
        found=False
        name=''.join(str(name).split()) #remove spaces in source names
        for target_name in target_names:
            target_name=''.join(target_name.split()) #remove spaces in target names just to be sure
            if target_name in name:
                found=True
                output=np.append(output,[name,df["UTC Obs Startyyyy-mm-dd hh:mm:ss"][ind],df["UTC Obs Endyyyy-mm-dd hh:mm:ss"][ind],df["PI"][ind],df["RAhh:mm:ss"][ind],df["DECdd:mm:ss"][ind]])
        
        #also check for FRB observations
        if "FRB" in name:
            output=np.append(output,[name,df["UTC Obs Startyyyy-mm-dd hh:mm:ss"][ind],df["UTC Obs Endyyyy-mm-dd hh:mm:ss"][ind],df["PI"][ind],df["RAhh:mm:ss"][ind],df["DECdd:mm:ss"][ind]])
        
        #check redshift of all scheduled sources to not miss ToO sources 
        try:
            average_redshift=np.average(Ned.get_table(name,table="redshifts")["Published Redshift"])
        except:
            average_redshift=1
        if average_redshift <0.2 and not found:
            output=np.append(output,["ToO " + name+" (z="+str(round(average_redshift,2))+")",df["UTC Obs Startyyyy-mm-dd hh:mm:ss"][ind],df["UTC Obs Endyyyy-mm-dd hh:mm:ss"][ind],df["PI"][ind],df["RAhh:mm:ss"][ind],df["DECdd:mm:ss"][ind]])
            
        
  
    return output

def sendEmail(password="IceCube!3",message="Testmail"):
    port = 587
    smtp_server="smtp.mail.de"
    sender_email="neutrino.alert@mail.de"
    
    context=ssl.create_default_context()
    
    with smtplib.SMTP(smtp_server, port) as server:
        server.starttls(context=context)
        server.login("neutrino.alert@mail.de",password)
        print("sending email")
        server.sendmail("neutrino.alert@mail.de","florian.eppel@gmx.de",message)

while True:
    data=checkObs()
    #generate message to send
    rows=[]
    for j in range(int(len(data)/6)):
        row=""
        for i in range(4):
            row=row+data[i+j*6]+" "
        obs_start=Time(data[1+j*6],format="iso",scale="utc")
        obs_end=Time(data[2+j*6],format="iso",scale="utc")
        ra=data[4+j*6]
        dec=data[5+j*6]

        coords=SkyCoord(ra+" "+dec,frame="icrs",unit=(u.hourangle,u.deg))
        
        exposure=(obs_end.mjd-obs_start.mjd)*24
        delta_time=np.linspace(0,exposure,100)*u.hour

        for telescope in telescopes:
            frame=AltAz(obstime=obs_start+delta_time,location=EarthLocation(lat=telescope["lat"],lon=telescope["lon"]))
        
            alt=coords.transform_to(frame).alt
            obs=obs_start+delta_time[np.where(alt>telescope["elev_lim"])[0]]
            if len(obs)>0:
                info=[np.min(obs),np.max(obs)]
                row=row+"\n "+telescope["name"] + ": " + info[0].iso + "-" + info[1].iso
        rows=np.append(rows,row)
        
        """

        frame_eff=AltAz(obstime=obs_start+delta_time,location=effelsberg)
        frame_parkes=AltAz(obstime=obs_start+delta_time,location=parkes)
        frame_fast=AltAz(obstime=obs_start+delta_time,location=fast)

        alt_eff=coords.transform_to(frame_eff).alt
        alt_parkes=coords.transform_to(frame_parkes).alt
        alt_fast=coords.transform_to(frame_fast).alt

        eff_obs=obs_start+delta_time[np.where(alt_eff>eff_elev_lim)[0]]
        parkes_obs=obs_start+delta_time[np.where(alt_parkes>parkes_elev_lim)[0]]
        fast_obs=obs_start+delta_time[np.where(alt_fast>fast_elev_lim)[0]]
        
        if len(eff_obs)>0:
            eff_info=[np.min(eff_obs),np.max(eff_obs)]
            row=row+"\n Effelsberg: "+eff_info[0].iso+"-"+eff_info[1].iso
        if len(parkes_obs)>0:
            parkes_info=[np.min(parkes_obs),np.max(parkes_obs)]
            row=row+"\n Parkes:"+parkes_info[0].iso+"-"+parkes_info[1].iso
        if len(fast_obs)>0:
            fast_info=[np.min(fast_obs),np.max(fast_obs)]
            row=row+"\n FAST:"+fast_info[0].iso+"-"+fast_info[1].iso
        
        rows=np.append(rows,row)
        
        """

    output_string=""
    for row in rows:
        output_string=output_string+row+"\n"+"\n"
    
    
    output_string=" Hi Flo,\n"+"Here is a list of upcoming XMM-Newton Targets to shadow search:\n\n"+output_string
    
    sendEmail(message=output_string)
    
    time.sleep(86400) #sleep for one day
