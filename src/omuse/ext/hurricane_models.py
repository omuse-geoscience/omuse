import datetime
import itertools

import numpy

from amuse.units import trigo
from omuse.units import units
from omuse.units import constants
from omuse.ext import spherical_geometry

try:
    from itertools import zip_longest
except ImportError:
    from itertools import izip_longest as zip_longest

# Best Track/Objective Aid/Wind Radii Format 
# http://www.nrlmry.navy.mil/atcf_web/docs/database/new/abrdeck.html

fields=["BASIN","CY","YYYYMMDDHH","TECHNUM/MIN","TECH","TAU","LatN/S","LonE/W", 
        "VMAX","MSLP","TY","RAD","WINDCODE","RAD1","RAD2","RAD3","RAD4","RADP",
        "RRP","MRD","GUSTS","EYE","SUBREGION","MAXSEAS","INITIALS","DIR","SPEED",
        "STORMNAME","DEPTH","SEAS","SEASCODE","SEAS1","SEAS2","SEAS3","SEAS4",
        "USERDEFINED","userdata"]

fields_description={
"BASIN" : "basin, e.g. WP, IO, SH, CP, EP, AL, SL",
"CY" : "annual cyclone number: 1 through 99",
"YYYYMMDDHH" : "Warning Date-Time-Group: 0000010100 through 9999123123. (note, 4 digit year)",
"TECHNUM/MIN" : "objective technique sorting number, minutes for best track: 00 - 99",
"TECH" : "acronym for each objective technique or CARQ or WRNG, BEST for best track.",
"TAU" : "forecast period: -24 through 240 hours, 0 for best-track, negative taus used for CARQ and WRNG records.",
"LatN/S" : "Latitude (tenths of degrees) for the DTG: 0 through 900, N/S is the hemispheric index.",
"LonE/W" : "Longitude (tenths of degrees) for the DTG: 0 through 1800, E/W is the hemispheric index.",
"VMAX" : "Maximum sustained wind speed in knots: 0 through 300.",
"MSLP" : "Minimum sea level pressure, 1 through 1100 MB.",
"TY" : """Level of tc development: 
DB - disturbance,
TD - tropical depression,
TS - tropical storm,
TY - typhoon,
ST - super typhoon,
TC - tropical cyclone,
HU - hurricane,
SD - subtropical depression,
SS - subtropical storm,
EX - extratropical systems,
IN - inland,
DS - dissipating,
LO - low,
WV - tropical wave,
ET - extrapolated,
XX - unknown. """,
"RAD" : "Wind intensity (kts) for the radii defined in this record: 34, 50, 64.",
"WINDCODE" : """Radius code:
AAA - full circle
QQQ - quadrant (NNQ, NEQ, EEQ, SEQ, SSQ, SWQ, WWQ, NWQ)""",
"RAD1" : "If full circle, radius of specified wind intensity, If semicircle or quadrant, radius of specified wind intensity of circle portion specified in radius code. 0 - 1200 nm.",
"RAD2" : "If full circle this field not used, If semicicle, radius (nm) of specified wind intensity for semicircle not specified in radius code, If quadrant, radius (nm) of specified wind intensity for 2nd quadrant (counting clockwise from quadrant specified in radius code). 0 through 1200 nm.",
"RAD3" : "If full circle or semicircle this field not used, If quadrant, radius (nm) of specified wind intensity for 3rd quadrant (counting clockwise from quadrant specified in radius code). 0 through 1200 nm.",
"RAD4" : "If full circle or semicircle this field not used, If quadrant, radius (nm) of specified wind intensity for 4th quadrant (counting clockwise from quadrant specified in radius code). 0 through 1200 nm.",
"RADP" : "pressure in millibars of the last closed isobar, 900 - 1050 mb.",
"RRP" : "radius of the last closed isobar in nm, 0 - 9999 nm.",
"MRD" : "radius of max winds, 0 - 999 nm.",
"GUSTS" : "gusts, 0 through 995 kts.",
"EYE" : "eye diameter, 0 through 999 nm.",
"SUBREGION" : """subregion code: W, A, B, S, P, C, E, L, Q.
A - Arabian Sea
B - Bay of Bengal
C - Central Pacific
E - Eastern Pacific
L - Atlantic
P - South Pacific (135E - 120W)
Q - South Atlantic
S - South IO (20E - 135E)
W - Western Pacific""",
"MAXSEAS" : "max seas: 0 through 999 ft.",
"INITIALS" : "Forecaster's initials, used for tau 0 WRNG, up to 3 chars.",
"DIR" : "storm direction in compass coordinates, 0 - 359 degrees.",
"SPEED" : "storm speed, 0 - 999 kts.",
"STORMNAME" : """literal storm name, NONAME or INVEST. TCcyx used pre-1999, where:
cy = Annual cyclone number 01 through 99
x = Subregion code: W, A, B, S, P, C, E, L, Q.
A - Arabain Sea
B - Bay of Bengal
C - Central Pacific
E - Eastern Pacific
L - Atlantic
P - South Pacific (135E - 120W)
Q - South Atlantic
S - South IO (20E - 135E)
W - Western Pacific""",
"DEPTH" : "system depth, D-deep, M-medium, S-shallow, X-unknown",
"SEAS" : "Wave height for radii defined in SEAS1-SEAS4, 0-99 ft.",
"SEASCODE" : """Radius code: 
AAA - full circle
QQQ - quadrant (NNQ, NEQ, EEQ, SEQ, SSQ, SWQ, WWQ, NWQ)""",
"SEAS1" : "first quadrant seas radius as defined by SEASCODE, 0 through 999 nm.",
"SEAS2" : "second quadrant seas radius as defined by SEASCODE, 0 through 999 nm.",
"SEAS3" : "third quadrant seas radius as defined by SEASCODE, 0 through 999 nm.",
"SEAS4" : "fourth quadrant seas radius as defined by SEASCODE, 0 through 999 nm.",
"USERDEFINED" : "20 character description of format to follow.",
"userdata" : "user data section as indicated by USERDEFINED parameter."
}

class ATCF_format_reader(object):
    """
    class to extract information from file in the
    
    Best Track/Objective Aid/Wind Radii Format 
    http://www.nrlmry.navy.mil/atcf_web/docs/database/new/abrdeck.html

    r=ATCF_format_reader(filename)
    """
    def __init__(self,filename="fort.22"):
        self.filename=filename
        self.data=None
        
    def read(self):
        """
        read data (store on data attribute)
        """
        f=open(self.filename,"r")
        lines=f.readlines()
        f.close()
        
        data=[]
        for line in lines:
          l=line[:-1].split(',')
          item=dict()
          for x,y in zip_longest(fields,l[:len(fields)], fillvalue=""):
              item[x]=y
          data.append(item)
        self.data=data
    def parse_datetime(self, s):
      return datetime.datetime.strptime(s.strip(),"%Y%m%d%H")
    def parse_lon(self,s):
      return ((s[-1]=="E")-(s[-1]=="W"))*float(s[:-1])/10 | units.deg
    def parse_lat(self,s):
      return ((s[-1]=="N")-(s[-1]=="S"))*float(s[:-1])/10 | units.deg
    def parse_rrp(self,s):
      return (0. if s=="" else float(s))  | units.nautical_mile
    def parse_mrd(self,s):
      return (0. if s=="" else float(s))  | units.nautical_mile
    def parse_vmax(self,s):
      return (0. if s=="" else float(s))  | units.knot
    def parse_mslp(self,s):
      return (0. if s=="" else float(s))  | units.mbar
    def parse_data(self):
      """
      parse data and return storm track
      (list of dicts with properties at different times)
      """
      last=""
      data=[]
      for dat in self.data:
        if dat["YYYYMMDDHH"]==last:
          continue
        
        if dat["TECH"].strip() != "BEST":
            raise Exception("forecast type not yet supported: {0}".format(dat["TECH"]))
        
        t=self.parse_datetime(dat["YYYYMMDDHH"])
        lon=self.parse_lon(dat["LonE/W"])
        lat=self.parse_lat(dat["LatN/S"])
        rrp=self.parse_rrp(dat["RRP"])
        mrd=self.parse_mrd(dat["MRD"])
        vmax=self.parse_vmax(dat["VMAX"])
        mslp=self.parse_mslp(dat["MSLP"])  
        
        if rrp.number == 0.:
          rrp=data[-1]["rrp"]
        if mrd.number == 0.:
          mrd=data[-1]["mrd"]
        if vmax.number == 0.:
          vmax=data[-1]["vmax"]
        if mslp.number == 0.:
          mslp=data[-1]["mslp"]
        
        data.append(dict(time=t,lon=lon,lat=lat,rrp=rrp,vmax=vmax,mslp=mslp,mrd=mrd))

        last=dat["YYYYMMDDHH"]
      return data

class HollandHurricane(object):
    """
    this class implements the Holland 1980 analytic wind and pressure
    profile. It can be initialized from a file in ATCF format or from 
    a track, which is to be provided as list of dicts with the 
    time, lon, lat, rrp (radius to last closed pressure contour), 
    vmax (maximum wind speed), mslp (maximum sea level pressure), 
    mrd (radius to max wind), example:
    
    h=HollandHurricane(nodes, "fort.22")
    """
    u10_boundary_layer_conversion_factor=0.7
    one_to_ten_minute_average_factor=0.88
    ambient_pressure=1013. | units.mbar
    air_density=1.15 | units.kg/units.m**3
    
    def __init__(self,nodes,file_or_track="fort.22"):
        if nodes:
            self.nodes=nodes.copy()
        else:
            self.nodes=None
        if isinstance(file_or_track,str):
            self.initialize_from_file(file_or_track)
        elif isinstance(file_or_track, list):
            self.initialize_from_track(file_or_track)
        else:
            raise Exception("provide filename or storm track")
        self.evolve_model(0. | units.s)
    def initialize_from_file(self,filename="fort.22"):
        r=ATCF_format_reader(filename=filename)
        r.read()
        track=r.parse_data()
        self.initialize_from_track(track)
    def initialize_from_track(self,track):
        self.add_translation_velocity(track)
        self.track=track
        t0=track[0]["time"]
        self.ndata=len(track)
        self.times=[(x["time"]-t0).total_seconds() for x in track] | units.s
    def add_translation_velocity(self,track):
        n=len(track)
        for i in range(n):
            i0=max(0,i-1)
            i1=min(n-1,i+1)
            t0=track[i0]["time"]
            lat0=track[i0]["lat"]
            lon0=track[i0]["lon"]
            t1=track[i1]["time"]
            lat1=track[i1]["lat"]
            lon1=track[i1]["lon"]
            dt=(t1-t0).total_seconds() | units.s
            d=spherical_geometry.distance(lat0,lon0,lat1,lon1, (1 | units.Rearth))
            v=d/dt
            dir_angle=trigo.arctan2(trigo.in_rad(lat1-lat0),trigo.in_rad(lon1-lon0))
            vx=v*trigo.cos(dir_angle)
            vy=v*trigo.sin(dir_angle)
            track[i]["vx"]=vx.in_(units.km/units.hour)
            track[i]["vy"]=vy.in_(units.km/units.hour)
    def interpolate_track(self,tend):
        index=numpy.interp(tend.value_in(units.s),self.times.value_in(units.s), numpy.arange(self.ndata))
        i=int(numpy.floor(index))
        dx=index-i
        i1=i+1
        if i==self.ndata-1:
          i1=self.ndata-1
        lat=(1-dx)*self.track[i]["lat"]+dx*self.track[i1]["lat"]
        lon=(1-dx)*self.track[i]["lon"]+dx*self.track[i1]["lon"]
        mslp=(1-dx)*self.track[i]["mslp"]+dx*self.track[i1]["mslp"]
        vmax=(1-dx)*self.track[i]["vmax"]+dx*self.track[i1]["vmax"]
        mrd=(1-dx)*self.track[i]["mrd"]+dx*self.track[i1]["mrd"]
        rrp=(1-dx)*self.track[i]["rrp"]+dx*self.track[i1]["rrp"]
        vx=(1-dx)*self.track[i]["vx"]+dx*self.track[i1]["vx"]
        vy=(1-dx)*self.track[i]["vy"]+dx*self.track[i1]["vy"]
        return lat,lon,mslp,vmax,mrd,rrp,vx,vy
    def evolve_model(self,tend):
        if not self.nodes:
            return
        lat,lon,mslp,vmax,mrd,rrp,vx,vy=self.interpolate_track(tend)
        
        central_pressure_deficit=self.ambient_pressure-mslp
        if central_pressure_deficit < 100| units.Pa:
            central_pressure_deficit=100. | units.Pa
        
        v_trans=(vx**2+vy**2)**0.5
        
        vmax=vmax-v_trans
        vmax=vmax/self.u10_boundary_layer_conversion_factor
        
        B=self.air_density*numpy.e*vmax**2/central_pressure_deficit
        if B<1: B=1
        if B>2.5: B=2.5
        
        coriolis_f=constants.coriolis_frequency(lat)
        r=spherical_geometry.distance(lat,lon,self.nodes.lat,self.nodes.lon, R=1.| units.Rearth)
        
        a=numpy.where(r.number==0.)
        r[a]=mrd*constants.eps
        
        theta=trigo.arctan2(trigo.in_rad(self.nodes.lat-lat),trigo.in_rad(self.nodes.lon-lon))
        self.nodes.pressure=mslp+central_pressure_deficit*numpy.exp(-(mrd/r)**B)
        
        v=((mrd/r)**B*numpy.exp(1-(mrd/r)**B)*vmax**2+r**2*coriolis_f**2/4)**0.5-r*coriolis_f/2
        
        v_trans_x=v/vmax*vx
        v_trans_y=v/vmax*vy
        
        fac=self.u10_boundary_layer_conversion_factor*self.one_to_ten_minute_average_factor
        
        self.nodes.vx=-v*trigo.sin(theta)*fac+v_trans_x
        self.nodes.vy=v*trigo.cos(theta)*fac+v_trans_y
