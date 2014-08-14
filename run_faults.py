#!/usr/bin/env python

import os

import batch

# import clawpack.geoclaw.okada2 as okada

class FaultJob(batch.Job):

    def __init__(self, deformation_file, topo_type=None):

        super(FaultJob, self).__init__()

        self.type = "compsyn"
        self.name = "guerrero_gap"
        self.executable = 'xgeoclaw'

        # Given an xyzt file, use this to run simulation
        self.deformation_file = os.path.abspath(deformation_file)
        self.prefix = os.path.basename(deformation_file).split('.')[0]

        # Data objects
        import setrun
        self.rundata = setrun.setrun()

        # Add deformation file
        self.rundata.dtopo_data.dtopofiles = []
        if topo_type is None:
            # Try to look at suffix for type
            extension = os.path.splitext(deformation_file)[1][1:]
            if extension[:2] == "tt":
                topo_type = int(extension[2])
            elif extension == 'xyz':
                topo_type = 1
            else:
                # Default to 3
                topo_type = 3
        self.rundata.dtopo_data.dtopofiles.append(
                                          [topo_type,3,5,self.deformation_file])
        # Add earthquake to regions list, force refinement for 1 minute
        regions = self.rundata.regiondata.regions
        regions.append([7, 7, 0.0, 60e0, -99.463937141374046, 
                                         -98.485815563597853, 
                                          16.622495995962375, 
                                          17.490586994378546 ])


    def __str__(self):
        output = super(FaultJob, self).__str__()
        if self.deformation_file == "okada.xyzt":
            output += "\n  Fault parameters: %s" % str(self.fault_params)
        else:
            output += "\n  Deformation File: %s" % self.deformation_file
        return output


if __name__ == '__main__':

    deformation_files = ["./bathy/rot_gapSzvT.xyzt",
                         "./bathy/rot_gapSvr1zvT.xyzt", 
                         "./bathy/rot_gapSvr2zvT.xyzt",
                         "./bathy/rot_west_gapSzvT.xyzt",
                         "./bathy/rot_west_gapSvr1zvT.xyzt", 
                         "./bathy/rot_west_gapSvr2zvT.xyzt",
                         "./bathy/okada_1957_du165.tt3",
                         "./bathy/okada_1957_du200.tt3"]

    tests = []
    for deformation_file in deformation_files:
        tests.append(FaultJob(deformation_file))
    
    controller = batch.BatchController(tests)
    print controller
    controller.tar = True
    controller.run()
