from datetime import datetime, tzinfo, timezone
from .momenttensor import moment_scale_factor
from .specfile import AMP_STYLE_DISP, AMP_STYLE_VEL

#import .earthmodel
DIST_SINGLE=1
DIST_REGULAR=0
DIST_IRREGULAR=-1

emptyStationXML = """
<?xml version="1.0" ?>
<FDSNStationXML xmlns="http://www.fdsn.org/xml/station/1" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.fdsn.org/xml/station/1 http://www.fdsn.org/xml/station/fdsn-station-1.1.xsd" schemaVersion="1.1">
    <Source>mgenkennett</Source>
    <Sender>pyreflect</Sender>
    <Created>{now}</Created>
    <Network code="{netcode}" startDate="{start}" >
        <Station code="{stacode}" startDate="{start}">
            <Latitude>{lat}</Latitude>
            <Longitude>{lon}</Longitude>
            <Elevation>0.0</Elevation>
            <Site>
               <Name>{sitename}</Name>
            </Site>
            <Channel code="{bandcode}{gaincode}Z" locationCode="{loccode}" startDate="{start}">
                <Latitude>{lat}</Latitude>
                <Longitude>{lon}</Longitude>
                <Elevation>0</Elevation>
                <Depth>0</Depth>
                <Azimuth>0</Azimuth>
                <Dip>-90</Dip>
                <Type>CONTINUOUS</Type>
                <Type>GEOPHYSICAL</Type>
                <SampleRate>{sps}</SampleRate>
                <Response>
                    <InstrumentSensitivity>
                        <Value>{gain}</Value>
                        <Frequency>0.0</Frequency>
                        <InputUnits>
                            <Name>{inputunits}</Name>
                        </InputUnits>
                        <OutputUnits>
                            <Name>counts</Name>
                            <Description>Digital Counts</Description>
                        </OutputUnits>
                    </InstrumentSensitivity>
                </Response>
            </Channel>
            <Channel code="{bandcode}{gaincode}R" locationCode="{loccode}" startDate="{start}">
                <Latitude>{lat}</Latitude>
                <Longitude>{lon}</Longitude>
                <Elevation>0</Elevation>
                <Depth>0</Depth>
                <Azimuth>{radialaz}</Azimuth>
                <Dip>0</Dip>
                <Type>CONTINUOUS</Type>
                <Type>GEOPHYSICAL</Type>
                <SampleRate>{sps}</SampleRate>
                <Response>
                    <InstrumentSensitivity>
                        <Value>{gain}</Value>
                        <Frequency>0.0</Frequency>
                        <InputUnits>
                            <Name>{inputunits}</Name>
                        </InputUnits>
                        <OutputUnits>
                            <Name>counts</Name>
                            <Description>Digital Counts</Description>
                        </OutputUnits>
                    </InstrumentSensitivity>
                </Response>
            </Channel>
            <Channel code="{bandcode}{gaincode}T" locationCode="{loccode}" startDate="{start}">
                <Latitude>{lat}</Latitude>
                <Longitude>{lon}</Longitude>
                <Elevation>0</Elevation>
                <Depth>0</Depth>
                <Azimuth>{transverseaz}</Azimuth>
                <Dip>0</Dip>
                <Type>CONTINUOUS</Type>
                <Type>GEOPHYSICAL</Type>
                <SampleRate>{sps}</SampleRate>
                <Response>
                    <InstrumentSensitivity>
                        <Value>{gain}</Value>
                        <Frequency>0.0</Frequency>
                        <InputUnits>
                            <Name>{inputunits}</Name>
                        </InputUnits>
                        <OutputUnits>
                            <Name>counts</Name>
                            <Description>Digital Counts</Description>
                        </OutputUnits>
                    </InstrumentSensitivity>
                </Response>
            </Channel>
        </Station>
    </Network>
</FDSNStationXML>
"""

def create_metadata(network, station, loccode, bandcode, gaincode, model, scalar_moment_N_m, ampStyle):
    inputunits = "m/s"
    if ampStyle == AMP_STYLE_DISP:
        inputunits = "m"
    data = {
      "netcode": network.code,
      "stacode": station.code,
      "loccode": loccode,
      "bandcode": bandcode,
      "gaincode": gaincode,
      "now": datetime.utcnow(),
      "start": station.start_date.format_iris_web_service(),
      "lat": station.latitude,
      "lon": station.longitude,
      "sitename": station.site.name,
      "radialaz": model.distance['azimuth'],
      "transverseaz": (model.distance['azimuth']+90) % 360,
      "gain": 1/moment_scale_factor(scalar_moment_N_m),
      "inputunits": inputunits,
      "sps": model.frequency['nyquist']
    }
    return emptyStationXML.format(**data).strip()


def create_fake_metadata(model, loccode, bandcode, gaincode, ampStyle=AMP_STYLE_VEL):
    network_code="XX"
    gain=1.0
    inputunits = "m/s"
    if ampStyle == AMP_STYLE_DISP:
        inputunits = "m"
    distList = model.list_distances()
    fakeXml = None
    for d in distList:
        deg = d/111.19
        station_code = str(round(deg, 2)).replace('.','_')
        data = {
          "netcode": network_code,
          "stacode": station_code,
          "loccode": loccode,
          "bandcode": bandcode,
          "gaincode": gaincode,
          "now": datetime.utcnow(),
          "start": "1900-01-01T00:00:00",
          "lat": 0.0,
          "lon": deg,
          "sitename": f"fake {deg} deg",
          "radialaz": model.distance['azimuth'],
          "transverseaz": (model.distance['azimuth']+90) % 360,
          "gain": gain,
          "inputunits": inputunits,
          "sps": model.frequency['nyquist']
        }
        if fakeXml is None:
            fakeXml = emptyStationXML.format(**data).strip()
        else:
            lines = fakeXml.split('\n')
            lines = lines[:-2]
            currlines = emptyStationXML.format(**data).strip().split('\n')
            currlines = currlines[6:]
            lines.extend(currlines)
            fakeXml = "\n".join(lines)
    return fakeXml
