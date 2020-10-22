import time
from pymongo import MongoClient
import gridfs
import base64
import io
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from datetime import datetime, timedelta, timezone
from siphon.simplewebservice.wyoming import WyomingUpperAir
from metpy.units import units
from metpy.plots import add_timestamp, SkewT, Hodograph
import metpy.calc as mpcalc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import scipy.ndimage as ndimage
from siphon.catalog import TDSCatalog
from matplotlib.colors import LinearSegmentedColormap

mpl.use("Agg")


client = MongoClient(os.environ["MONGODB_CLIENT"])
db = client.wx_gfx
fs = gridfs.GridFS(db)


def plot_skewt(df):
    hght = df["height"].values * units.hPa
    p = df["pressure"].values * units.hPa
    T = df["temperature"].values * units.degC
    Td = df["dewpoint"].values * units.degC
    wind_speed = df["speed"].values * units.knots
    wind_dir = df["direction"].values * units.degrees
    u, v = mpcalc.wind_components(wind_speed, wind_dir)

    fig = plt.figure(figsize=(9, 12))
    skew = SkewT(fig, rotation=45)

    skew.plot(p, T, "r")
    skew.plot(p, Td, "g")
    skew.plot_barbs(p, u, v)
    skew.ax.set_ylim(1000, 100)
    skew.ax.set_xlim(-40, 60)

    lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
    skew.plot(lcl_pressure, lcl_temperature, "ko", markerfacecolor="black")

    prof = mpcalc.parcel_profile(p, T[0], Td[0]).to("degC")
    skew.plot(p, prof, "k", linewidth=2)

    skew.ax.axvline(0, color="c", linestyle="--", linewidth=2)

    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()

    ax_hod = inset_axes(skew.ax, "40%", "40%", loc=2)
    h = Hodograph(ax_hod, component_range=80.0)
    h.add_grid(increment=20)
    h.plot_colormapped(u, v, hght)

    return skew


def make_name(site, date, time):
    if date:
        return "{site}_{dt:%Y%m%d_%H%M}.png".format(site=site, dt=time)
    else:
        return "{site}.svg".format(site=site)


def generate_sounding_plot(site, date=None):

    if date:
        request_time = datetime.strptime(date, "%Y%m%d%H")
    else:
        now = datetime.now(timezone.utc) - timedelta(hours=2)
        request_time = now.replace(
            hour=(now.hour // 12) * 12, minute=0, second=0
        )

    df = WyomingUpperAir.request_data(request_time, site)
    skewt = plot_skewt(df)

    add_timestamp(
        skewt.ax, request_time, y=1.02, x=0, ha="left", fontsize="large"
    )
    skewt.ax.set_title(site)

    bio = io.BytesIO()
    skewt.ax.figure.savefig(bio, format="svg")
    bio.seek(0)
    b64 = base64.b64encode(bio.read())
    try:
        file = fs.find_one({"filename": site})
        fs.delete(file._id)
    except Exception:
        pass
    fs.put(b64, filename=site, timestamp=datetime.now(timezone.utc))

    plt.close(skewt)


station_list = [
    "OAK",
    "REV",
    "LKN",
    "SLC",
    "GJT",
    "DNR",
    "VBG",
    "EDW",
    "DRA",
    "FGZ",
    "ABQ",
    "AMA",
    "NKX",
    "TUS",
    "EPZ",
    "MAF",
    "FWD",
    "SHV",
    "DRT",
    "CRP",
    "BRO",
    "LCH",
    "SHV",
    "OUN",
    "DDC",
    "RIW",
    "BOI",
    "MFR",
    "SLE",
    "UIL",
    "OTX",
    "TFX",
    "GGW",
    "UNR",
    "BIS",
    "ABR",
    "OAX",
    "TOP",
    "SGF",
    "LZK",
    "SIL",
    "INL",
    "MPX",
    "DVN",
    "ILX",
    "JAN",
    "APX",
    "DTX",
    "ILN",
    "BNA",
    "BMX",
    "FFC",
    "TLH",
    "PIT",
    "RNK",
    "GSO",
    "CHS",
    "JAX",
    "TBW",
    "XMR",
    "MFL",
    "EYW",
    "MHX",
    "IAD",
    "BUF",
    "WAL",
    "ALB",
    "OKX",
    "GYX",
    "CHH",
    "CAR",
]


def generate_rap_plots():

    colors = ["#cbff30", "#7fff30", "#30ff56", "#129e10"]
    cm = LinearSegmentedColormap.from_list("cm", colors, N=len(colors) * 10)

    now = datetime.now(timezone.utc)

    # Define datasets to wanted.
    prop = "Relative_humidity_isobaric"

    # Use TDSCatalog to begin data access, NCSS to subset, & grab latest file.
    rap_cat = TDSCatalog(
        "https://thredds.unidata.ucar.edu/thredds/catalog/\
        grib/NCEP/RAP/CONUS_13km/latest.xml"
    )
    latestrap = rap_cat.datasets[0]
    ncss = latestrap.subset()

    # Query prop data.
    prop_data = ncss.query()
    prop_data.variables(prop)
    prop_data.add_lonlat().lonlat_box(north=59, south=15, east=-57, west=-139)
    prop_data.time(now)
    prop_dataq = ncss.get_data(prop_data)

    propc = prop_dataq.variables[prop]

    levels = dict(_925=33, _850=28, _700=24, _500=16, _300=8, _250=6, _200=4)

    for level in levels:
        name = level[1:]
        slice = levels[level]

        propcs = propc[0, slice, :, :]

        fnl_prop = ndimage.gaussian_filter(propcs, sigma=0.5, order=0)

        # Extract the lon/lat,
        lon = prop_dataq.variables["lon"][:]
        lat = prop_dataq.variables["lat"][:]

        # Create figure.
        fig = plt.figure(figsize=(10, 15))

        # Define projection.
        lc = ccrs.LambertConformal(
            central_longitude=-97.5,
            central_latitude=35,
            standard_parallels=(30, 60),
        )
        ax = fig.add_subplot(1, 1, 1, projection=lc)
        ax.set_extent([-123, -73, 23, 52], crs=ccrs.PlateCarree())

        # Create the map.
        ax.add_feature(
            cfeature.OCEAN.with_scale("50m"),
            facecolor="#34cceb",
            edgecolor="none",
            zorder=5,
        )
        ax.add_feature(
            cfeature.LAND.with_scale("50m"),
            edgecolor="dimgray",
            facecolor="#ede6af",
            zorder=0,
        )
        ax.add_feature(cfeature.BORDERS.with_scale("50m"), zorder=8)
        ax.add_feature(
            cfeature.LAKES.with_scale("50m"),
            linewidth=0.5,
            facecolor="#34cceb",
            edgecolor="dimgray",
            zorder=3,
        )
        ax.add_feature(
            cfeature.STATES.with_scale("50m"),
            linewidth=0.5,
            edgecolor="black",
            zorder=8,
        )

        # Plot prop data.
        cntr_prop = np.arange(50, 101, 1)
        try:
            ax.contourf(
                lon,
                lat,
                fnl_prop,
                cntr_prop,
                cmap=cm,
                zorder=7,
                transform=ccrs.PlateCarree(),
            )
        except Exception:
            pass

        ax.gridlines(zorder=9)

        bio = io.BytesIO()
        ax.figure.savefig(bio, format="svg", bbox_inches="tight")
        bio.seek(0)
        b64 = base64.b64encode(bio.read())
        try:
            file = fs.find_one({"filename": name})
            fs.delete(file._id)
        except Exception:
            pass
        fs.put(b64, filename=name, timestamp=datetime.now(timezone.utc))

        fig.clear()
        plt.close(fig)


if __name__ == "__main__":
    next_hour = datetime.now(timezone.utc)
    while True:
        if datetime.now(timezone.utc) >= next_hour:
            for station in station_list:
                try:
                    generate_sounding_plot(station)
                except Exception:
                    pass
            generate_rap_plots()

            print("got metpy")
            next_hour = datetime.now(timezone.utc) + timedelta(hours=1)
        else:
            print("skipping updates")
        time.sleep(60 * 10)
