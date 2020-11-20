from datetime import datetime, tzinfo, timezone
from .momenttensor import moment_scale_factor
from .specfile import AMP_STYLE_DISP, AMP_STYLE_VEL

emptyStationXML = """
<?xml version="1.0" encoding="ISO-8859-1"?>
<FDSNStationXML xmlns="http://www.fdsn.org/xml/station/1" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.fdsn.org/xml/station/1 http://www.fdsn.org/xml/station/fdsn-station-1.1.xsd" schemaVersion="1.1">
    <Source>mgenkennett</Source>
    <Sender>pyreflect</Sender>
    <Created>{now}</Created>
    <Network code="{netcode}" startDate="{start}" >
        <Station code="{stacode}" startDate="{start}">
            <Latitude>36.9488</Latitude>
            <Longitude>-116.2495</Longitude>
            <Elevation>1600.0</Elevation>
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

def createMetadata(network, station, loccode, bandcode, gaincode, model, scalar_moment_N_m, ampStyle):
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
