#!/usr/bin/env python
__author__ = 'adamkoziol,mikeknowles'


class MetadataReader(object):

    def reader(self):
        import os
        import json
        from accessoryFunctions import GenObject, MetadataObject
        for sample in self.metadata:
            metadatafile = os.path.join(self.path, sample.name, sample.name + "_metadata.json")
            if os.path.isfile(metadatafile):
                with open(metadatafile) as metadatareport:
                    jsondata = json.load(metadatareport)
                # Create the metadata objects
                metadata = MetadataObject()
                # Set the name
                metadata.name = sample.name
                # Initialise the metadata categories as GenObjects created using the appropriate key
                for attr in jsondata:
                    if not isinstance(jsondata[attr], dict):
                        setattr(metadata, attr, jsondata[attr])
                    else:
                        setattr(metadata, attr, GenObject(jsondata[attr]))
                # metadata.run = GenObject(jsondata['run'])
                # metadata.general = GenObject(jsondata['general'])
                # metadata.commands = GenObject(jsondata['commands'])
                self.samples.append(metadata)

    def __init__(self, inputobject):
        # self.metadata = inputobject.runmetadata.samples
        self.metadata = inputobject.samples
        self.path = inputobject.path
        self.samples = []
        self.reader()
