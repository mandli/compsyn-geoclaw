#!/usr/bin/env python

import os
import sys
import glob

import batch

class FaultJob(batch.Job):

    def __init__(self, deformation_file):

        super(FaultJob, self).__init__()

        self.type = "compsyn"
        self.name = "guerrero_gap"
        self.prefix = os.path.basename(deformation_file).split('.')[0]
        self.deformation_file = os.path.abspath(deformation_file)
        self.executable = 'xgeoclaw'

        # Data objects
        import setrun
        self.rundata = setrun.setrun()

        # Add deformation file
        self.rundata.dtopo_data.dtopofiles = []
        self.rundata.dtopo_data.dtopofiles.append([1,3,5,self.deformation_file])


    def __str__(self):
        output = super(FaultJob, self).__str__()
        output += "\n  Deformation File: %s" % self.deformation_file
        return output



if __name__ == '__main__':

    deformation_files = ["./bathy/rot_gapSzvT.xyzt",
                         "./bathy/rot_gapSvr1zvT.xyzt", 
                         "./bathy/rot_gapSvr2zvT.xyzt",
                         "./bathy/rot_west_gapSzvT.xyzt",
                         "./bathy/rot_west_gapSvr1zvT.xyzt", 
                         "./bathy/rot_west_gapSvr2zvT.xyzt"]

    # deformation_files = ["./bathy/rot_gapSzvT.xyzt",
    #                      "./bathy/rot_west_gapSzvT.xyzt"]
    # if len(sys.argv) > 1:
    #     deformation_files = sys.argv[1:]
    # else:
    #     deformation_files = glob.glob('./bathy/rot_gapS*.xyzt')
    #     for west_file in glob.glob('./bathy/rot_west_gapS*.xyzt'):
    #         deformation_files.append(west_file)

    tests = []
    for deformation_file in deformation_files:
        tests.append(FaultJob(deformation_file))
    
    controller = batch.BatchController(tests)
    print controller
    controller.tar = True
    controller.run()
