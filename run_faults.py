#!/usr/bin/env python

import os
import sys
import glob

import clawpack.tests as clawtests

class FaultTest(clawtests.Test):

    def __init__(self, deformation_file):

        super(FaultTest, self).__init__()

        self.type = "compsys"
        self.name = "guerrero_gap"
        self.prefix = os.path.basename(deformation_file).split('.')[0]
        self.deformation_file = os.path.abspath(deformation_file)
        self.executable = 'xgeoclaw'

        # Data objects
        import setrun
        self.rundata = setrun.setrun()

        # Add deformation file
        self.rundata.dtopo_data.dtopofiles = []
        self.rundata.dtopo_data.dtopofiles.append([1,5,5,self.deformation_file])


    def __str__(self):
        output = super(FaultTest, self).__str__()
        output += "\n  Deformation File: %s" % self.deformation_file
        return output



if __name__ == '__main__':

    if len(sys.argv) > 1:
        deformation_files = sys.argv[1:]
    else:
        deformation_files = glob.glob('./bathy/rot_gap*z*.xyzt')

    tests = []
    for deformation_file in deformation_files:
        tests.append(FaultTest(deformation_file))
    
    controller = clawtests.TestController(tests)
    print controller
    controller.tar = True
    controller.run()
