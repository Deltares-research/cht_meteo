
import json
import numpy as np

def write_wind_to_json(dataset, file_name, time_range=None, iref=1, js=False):

    # Convert numpy.datetime64 to datetime
    time = dataset.ds.time.values.astype("M8[s]").astype("O")
    x    = dataset.ds.lon.values
    y    = dataset.ds.lat.values

    if not time_range:
        time_range = []
        time_range.append(time[0])
        time_range.append(time[-1])

    data = []

    header = {
        "discipline": 0,
        "disciplineName": "Meteorological products",
        "gribEdition": 2,
        "gribLength": 76420,
        "center": 7,
        "centerName": "US National Weather Service - NCEP(WMC)",
        "subcenter": 0,
        "refTime": "2016-04-30T06:00:00.000Z",
        "significanceOfRT": 1,
        "significanceOfRTName": "Start of forecast",
        "productStatus": 0,
        "productStatusName": "Operational products",
        "productType": 1,
        "productTypeName": "Forecast products",
        "productDefinitionTemplate": 0,
        "productDefinitionTemplateName": "Analysis/forecast at horizontal level/layer at a point in time",
        "parameterCategory": 2,
        "parameterCategoryName": "Momentum",
        "parameterNumber": 2,
        "parameterNumberName": "U-component_of_wind",
        "parameterUnit": "m.s-1",
        "genProcessType": 2,
        "genProcessTypeName": "Forecast",
        "forecastTime": 0,
        "surface1Type": 103,
        "surface1TypeName": "Specified height level above ground",
        "surface1Value": 10.0,
        "surface2Type": 255,
        "surface2TypeName": "Missing",
        "surface2Value": 0.0,
        "gridDefinitionTemplate": 0,
        "gridDefinitionTemplateName": "Latitude_Longitude",
        "numberPoints": 65160,
        "shape": 6,
        "shapeName": "Earth spherical with radius of 6,371,229.0 m",
        "gridUnits": "degrees",
        "resolution": 48,
        "winds": "true",
        "scanMode": 0,
        "nx": 360,
        "ny": 181,
        "basicAngle": 0,
        "subDivisions": 0,
        "lo1": 0.0,
        "la1": 90.0,
        "lo2": 359.0,
        "la2": -90.0,
        "dx": 1.0,
        "dy": 1.0,
    }

    header["lo1"] = float(min(x) + 360.0)
    header["lo2"] = float(max(x) + 360.0)
    header["la1"] = float(max(y))
    header["la2"] = float(min(y))
    header["dx"] = float(x[1] - x[0])
    header["dy"] = float(y[1] - y[0])
    header["nx"] = len(x)
    header["ny"] = len(y)
    header["numberPoints"] = len(x) * len(y)

    header_u = header.copy()
    header_v = header.copy()

    header_u["parameterNumberName"] = "U-component_of_wind"
    header_u["parameterNumber"] = 2
    header_v["parameterNumberName"] = "V-component_of_wind"
    header_v["parameterNumber"] = 3

    for it, t in enumerate(time):
        if t >= time_range[0] and t <= time_range[1]:
            dd = []

            tstr = t.strftime("%Y-%m-%dT%H:%M:%SZ")

            u_list = (
                np.flipud(np.around(dataset.ds["wind_u"].values[it, :, :], decimals=1))
                .flatten()
                .tolist()
            )
            data0 = {"header": header_u.copy(), "data": u_list}
            data0["header"]["refTime"] = tstr
            dd.append(data0)

            v_list = (
                np.flipud(np.around(dataset.ds["wind_v"].values[it, :, :], decimals=1))
                .flatten()
                .tolist()
            )
            data0 = {"header": header_v.copy(), "data": v_list}
            data0["header"]["refTime"] = tstr
            dd.append(data0)

            data.append(dd)

    json_string = json.dumps(data, separators=(",", ":"))
    fid = open(file_name, "w")
    if js:
        fid.write("wind = ")
    fid.write(json_string)
    fid.close()
