"""Generate a microscopic watertight tetrahedron and controls, using only Python's standard library."""
from pathlib import Path
import json
import zipfile
import xml.etree.ElementTree as ET

OUT = Path(__file__).resolve().parent
CORE = 'http://schemas.microsoft.com/3dmanufacturing/core/2015/02'
CONFIG = {'brim_type': 'no_brim', 'curr_bed_type': 'Textured PEI Plate', 'enable_prime_tower': '0', 'enable_support': '0', 'filament_settings_id': ['Generic PLA @BBL A1'], 'initial_layer_print_height': '0.2', 'layer_height': '0.2', 'nozzle_diameter': ['0.4'], 'print_settings_id': '0.20mm Standard @BBL A1', 'printable_height': '256', 'printer_settings_id': 'Bambu Lab A1 0.4 nozzle', 'raft_layers': '2', 'sparse_infill_pattern': 'grid'}

TRIANGLES = [(0,2,1),(0,1,3),(1,2,3),(2,0,3)]

def save(name, height):
    model = ET.Element('model', {'unit':'millimeter', 'xmlns':CORE})
    ET.SubElement(model, 'metadata', {'name':'Application'}).text = 'OrcaSlicer-2.4.2'
    ET.SubElement(model, 'metadata', {'name':'BambuStudio:3mfVersion'}).text = '1'
    mesh = ET.SubElement(ET.SubElement(ET.SubElement(model, 'resources'), 'object',
                        {'id':'1','type':'model','name':name}), 'mesh')
    width = 0.0002 if height < 0.01 else 10.
    vertices = [(0,0,0),(width,0,0),(0,width,0),(0,0,height)]
    vs = ET.SubElement(mesh, 'vertices')
    for x,y,z in vertices:
        ET.SubElement(vs,'vertex',{'x':str(x),'y':str(y),'z':str(z)})
    ts = ET.SubElement(mesh, 'triangles')
    for a,b,c in TRIANGLES:
        ET.SubElement(ts,'triangle',{'v1':str(a),'v2':str(b),'v3':str(c)})
    ET.SubElement(ET.SubElement(model,'build'),'item',
                  {'objectid':'1','transform':'1 0 0 0 1 0 0 0 1 123 123 0'})
    settings = ET.Element('config')
    obj = ET.SubElement(settings,'object',{'id':'1','auto_drop':'1'})
    ET.SubElement(obj,'metadata',{'key':'name','value':name})
    ET.SubElement(obj,'metadata',{'key':'extruder','value':'1'})
    part = ET.SubElement(obj,'part',{'id':'1','subtype':'normal_part'})
    ET.SubElement(part,'metadata',{'key':'name','value':name})
    plate = ET.SubElement(settings,'plate')
    ET.SubElement(plate,'metadata',{'key':'plater_id','value':'1'})
    instance = ET.SubElement(plate,'model_instance')
    for key,value in [('object_id','1'),('instance_id','0'),('identify_id','1')]:
        ET.SubElement(instance,'metadata',{'key':key,'value':value})
    config = CONFIG
    types = '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/><Default Extension="config" ContentType="text/xml"/></Types>'
    rels = '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Target="/3D/3dmodel.model" Id="rel0" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/></Relationships>'
    with zipfile.ZipFile(OUT/(name+'.3mf'),'w',zipfile.ZIP_DEFLATED) as z:
        z.writestr('[Content_Types].xml',types)
        z.writestr('_rels/.rels',rels)
        z.writestr('3D/3dmodel.model',ET.tostring(model,encoding='utf-8',xml_declaration=True))
        z.writestr('Metadata/model_settings.config',ET.tostring(settings,encoding='utf-8',xml_declaration=True))
        z.writestr('Metadata/project_settings.config',json.dumps(config))

if __name__ == '__main__':
    save('empty-layer-crash',0.002)
    save('normal-tetrahedron-control',1.)
    CONFIG['raft_layers'] = '0'
    save('no-raft-control',0.002)
