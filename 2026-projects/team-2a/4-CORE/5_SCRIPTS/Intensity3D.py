#! /usr/bin/python

import sys, os
import numpy
import scipy
import scipy.interpolate
import matplotlib
import matplotlib.pyplot as plt

USAGE = '''
USAGE: %s [Path to 3D dose file] [x|y|z] [depth]

The following example plots the z plane at z=10mm.

   python %s output.dat z 10.0
'''

class Intensity3DFileException(Exception):

    def __init__(self, message):
        self.message = message

class Intensity3D:

    ###################################################
    # PRIVATE CLASS METHODS
    ###################################################


    def get_num_bins(self, records):
        bin_lengths = records[0].split(",")

        if len(bin_lengths) != 3:
            raise Intensity3DFileException("Error: Incorrect bin lengths in Intensity3D file.")

        # python3: converted map to list
        num_bins = list(map(int, bin_lengths))

        return num_bins


    def get_bin_edges(self, records, num_bins):

        edges = [r.split(",") for r in records[1:4]]

        for i in [0,1,2]:
            if len(edges[i]) - 1 != num_bins[i]:
                err_msg = "Error:  wrong number of bin edges (%d) for %d bins in dimension %d" % (len(edges[i]), num_bins[i], i)
                print (err_msg)
                raise Intensity3DFileException(err_msg)

            # python3: converted map to list
            edges[i] = list(map(float, edges[i]))

        return edges


    def get_bin_values(self, records, num_bins):

        #get rid of trailing EOL characters
        while( len(records[-1].strip()) == 0):
            records.pop()

        if num_bins[2] != len(records[4:]):
            err_msg = "Error: number of records (%d) does not equal length of z dimension (%d)" % (num_bins[2], len(records[4:]))
            print (err_msg)
            raise Intensity3DFileException(err_msg)

        step_size = num_bins[0]*num_bins[1]

        array_3D = numpy.zeros(shape=(num_bins[0],num_bins[1],num_bins[2]))
        for k, rec in enumerate(records[4:]):
            rec = rec.strip()
            if len(rec) == 0:
                continue
            values = [float(value.strip()) for value in rec.split(",")]
            for j in range(num_bins[1]):
                for i in range(num_bins[0]):
                    bin_number = j*num_bins[0] + i
                    array_3D[i,j,k] = values[bin_number]

        if num_bins[0] == 1:
            array_3D_doubled = numpy.zeros(shape=(2, num_bins[1], num_bins[2]))
            array_3D_doubled[0,:,:] = array_3D[0,:,:]
            array_3D_doubled[1,:,:] = array_3D[0,:,:]
            array_3D = array_3D_doubled
        if num_bins[1] == 1:
            array_3D_doubled = numpy.zeros(shape=(num_bins[0], 2, num_bins[2]))
            array_3D_doubled[:,0,:] = array_3D[:,0,:]
            array_3D_doubled[:,1,:] = array_3D[:,0,:]
            array_3D = array_3D_doubled
        if num_bins[2] == 1:
            array_3D_doubled = numpy.zeros(shape=(num_bins[0],num_bins[1],2))
            array_3D_doubled[:,:,0] = array_3D[:,:,0]
            array_3D_doubled[:,:,1] = array_3D[:,:,0]
            array_3D = array_3D_doubled

        return array_3D


    def read_dose_file(self, file_path):

        #we need get rid of spaces and make commas the delimiter
        scrubber = lambda r: r.strip().strip(",").replace(", ", ",").replace(" ", ",")
        # python3: converted map to list
        dose_records = list(map(scrubber, open(file_path, 'r').readlines()))

        self.num_bins = self.get_num_bins(dose_records)
        self.edges = self.get_bin_edges(dose_records, self.num_bins)
        self.values = self.get_bin_values(dose_records, self.num_bins)

        return self.num_bins, self.edges, self.values


    def get_bin_centers(self, edges_array):

        lengths = [len(a) for a in edges_array]
        lprod = numpy.prod(lengths)
        number_of_vals = numpy.prod([len(a) -1 for a in edges_array])
        centers_array = []

        shape = []
        for i, edges in enumerate(edges_array):
            assert(min(edges) == edges[0])
            assert(max(edges) == edges[-1])

            #if there are only 2 edges then keep them both as the centers
            # because interpolations fails if any dimension has only 1
            if len(edges) == 2:
                centers_array.append(edges)
            else:
                step_size = (edges[-1] - edges[0])/float(len(edges) - 1)
                centers = numpy.array(edges[0:-1])
                centers +=  step_size * 0.5

                centers_array.append(centers)

        return centers_array

    def get_ranges(self):
        return [[min(self.bin_centers[0]), max(self.bin_centers[0])], \
                [min(self.bin_centers[1]), max(self.bin_centers[1])],
                [min(self.bin_centers[2]), max(self.bin_centers[2])]]

    def _validate_profile_coords(self, dimension, fixed_coordinates):
        """
        Validate the fixed coordinates used to build a 1D profile along `dimension`.
        Raises ValueError with a descriptive message if any fixed coordinate lies outside
        the available bin center ranges for that axis.
        """
        assert(len(fixed_coordinates) == 2)
        xyz = [0,1,2]
        xyz.remove(dimension)

        # build readable ranges
        ranges = self.get_ranges()
        out_of_bounds = []
        labels = ['x', 'y', 'z']
        for idx, coord in enumerate(fixed_coordinates):
            axis = xyz[idx]
            axis_min, axis_max = ranges[axis]
            if coord < axis_min or coord > axis_max:
                out_of_bounds.append((labels[axis], coord, axis_min, axis_max))

        if out_of_bounds:
            msgs = []
            for axis_label, coord, amin, amax in out_of_bounds:
                msgs.append(f"{axis_label} coordinate {coord} outside range [{amin}, {amax}]")
            raise ValueError("; ".join(msgs))

    def _validate_plane_coord(self, plane_index, depth):
        """
        Validate that the requested plane depth lies within the bin center range
        for the specified axis. Raises ValueError with descriptive message if not.
        """
        labels = ['x', 'y', 'z']
        amin = min(self.bin_centers[plane_index])
        amax = max(self.bin_centers[plane_index])
        if depth < amin or depth > amax:
            raise ValueError(
                f"Requested plane depth {depth} for axis {labels[plane_index]} "
                f"is outside bin-center range [{amin}, {amax}]. Check units/order."
            )

    #
    def getIntensityPlaneRange(self, plane_index, depth, range_min, range_max):

        assert(len(self.bin_centers) == 3)
        xyz = [0,1,2]
        xyz.remove(plane_index)
        labels = ['x', 'y', 'z']

        # Validate requested plane coordinate before building points
        self._validate_plane_coord(plane_index, depth)

        centers = list(self.bin_centers)
        centers[plane_index] = depth

        points = []
        h_vals = centers[xyz[0]]
        h_vals = h_vals[h_vals < range_max]
        h_vals = h_vals[h_vals > range_min]
        v_vals = centers[xyz[1]]
        v_vals = v_vals[v_vals < range_max]
        v_vals = v_vals[v_vals > range_min]
        for j in range(len(v_vals)):
            for i in range(len(h_vals)):
                point = [0,0,0]
                point[xyz[0]] = h_vals[i]
                point[xyz[1]] = v_vals[j]
                point[plane_index] = depth
                points.append(point)

        # defensive call to interpolator: catch and re-raise with helpful context
        try:
            values = self.f3D(points)
        except ValueError as e:
            pmin = numpy.min(points, axis=0)
            pmax = numpy.max(points, axis=0)
            grid_ranges = [(min(a), max(a)) for a in self.bin_centers]
            raise ValueError(
                f"Interpolator error: {e};\n"
                f"points bbox min={pmin}, max={pmax}; grid ranges={grid_ranges}"
            ) from e
        values = values.reshape(len(v_vals), len(h_vals))

        return values


    def getIntensityPlane(self, plane_index, depth):
        assert (len(self.bin_centers) == 3)
        xyz = [0, 1, 2]
        xyz.remove(plane_index)
        labels = ['x', 'y', 'z']

        # Validate the requested plane coordinate
        self._validate_plane_coord(plane_index, depth)

        centers = list(self.bin_centers)
        centers[plane_index] = depth

        points = []
        for j in range(len(centers[xyz[1]])):
            for i in range(len(centers[xyz[0]])):
                point = [0, 0, 0]
                point[xyz[0]] = centers[xyz[0]][i]
                point[xyz[1]] = centers[xyz[1]][j]
                point[plane_index] = depth
                points.append(point)

        # defensive call to interpolator with contextual error if it fails
        try:
            values = self.f3D(points)
        except ValueError as e:
            pmin = numpy.min(points, axis=0)
            pmax = numpy.max(points, axis=0)
            grid_ranges = [(min(a), max(a)) for a in self.bin_centers]
            raise ValueError(
                f"Interpolator error in getIntensityPlane: {e};\n"
                f"points bbox min={pmin}, max={pmax}; grid ranges={grid_ranges}"
            ) from e
        values = values.reshape(len(centers[xyz[1]]), len(centers[xyz[0]]))

        return values


    def get_profile(self, dimension, fixed_coordinates = [0.0, 0.0]):
        assert(len(fixed_coordinates) == 2)
        xyz = [0,1,2]
        xyz.remove(dimension)
        labels = ['x', 'y', 'z']

        centers = self.bin_centers[dimension]
        #ZZ = numpy.linspace(centers[0], centers[-1], num=5.0*len(centers), endpoint=True)
        ZZ = numpy.linspace(centers[0], centers[-1], num = int(5.0*len(centers)), endpoint=True)  # fixing error: TypeError: 'float' object cannot be interpreted as an integer

        # Validate fixed coordinates before building interpolation points
        self._validate_profile_coords(dimension, fixed_coordinates)

        points = []
        for p in ZZ:
            point = [0,0,0]
            point[xyz[0]] = fixed_coordinates[0]
            point[xyz[1]] = fixed_coordinates[1]
            point[dimension] = p
            points.append(point)

        values = self.f3D(points)

        return ZZ, values
    #
    def get_plot_range(self, plane_index, depth, range_min, range_max, title=" ", is_normalized=False):

        print ("get_plot_range", plane_index, depth, range_min, range_max)
        plane = self.getIntensityPlaneRange(plane_index, depth, range_min, range_max)
        plane /= numpy.max(plane.flatten())

        axes = [['x', self.edges[0]], ['y', self.edges[1]], ['z', self.edges[2]]]
        plane_string = axes[plane_index][0]
        del axes[plane_index]

        plt.clf()
        title = "%s = %.1f mm" % (plane_string, depth)
        plt.title(title, fontsize=22)
        plt.grid(True)
        # plt.xlabel("%s [mm]" % axes[0][0], fontsize=20)
        # plt.ylabel("%s [mm]" % axes[1][0], fontsize=20)

        # extent = [min(axes[0][1]), max(axes[0][1]), min(axes[1][1]), max(axes[1][1])]
        extent = [range_min, range_max, range_min, range_max]
        #plt.axis((min(axes[0][1]), max(axes[0][1]), min(axes[1][1]), max(axes[1][1])))
        plt.axis([range_min, range_max, range_min, range_max])

        plot = plt.imshow(plane[::-1, :], extent=extent)

        cb = plt.colorbar(plot)
        # cb.set_ticks(ticks = [0.0,0.25, 0.5, 0.75, 1.0, 1.25])

        filename = "%s/range_%s%.1f.png" \
                   % (self.output_folder, plane_string, depth)
        plt.gcf().set_size_inches(9, 6)
        plt.savefig(filename)

        return plot
        

    def get_plot(self, plane_index, depth, title=" ", is_normalized=False):
        plane = self.getIntensityPlane(plane_index, depth)

        if is_normalized:
            plane /= numpy.max(plane.flatten())

        axes = [['x', self.edges[0]], ['y', self.edges[1]], ['z', self.edges[2]]]
        plane_string = axes[plane_index][0]
        del axes[plane_index]

        plt.clf()
        title = "%s = %.1f mm" % (plane_string, depth)
        plt.title(title, fontsize=22)
        plt.grid(True)
        plt.xlabel("%s [mm]" % axes[0][0], fontsize=20)
        plt.ylabel("%s [mm]" % axes[1][0], fontsize=20)

        extent = [min(axes[0][1]), max(axes[0][1]), min(axes[1][1]), max(axes[1][1])]
        plt.axis(( min(axes[0][1]), max(axes[0][1]), min(axes[1][1]), max(axes[1][1]) ))

        plot = plt.imshow(plane[::-1,:], extent=extent)

        cb = plt.colorbar(plot)
        # cb.set_ticks(ticks = [0.0,0.25, 0.5, 0.75, 1.0, 1.25])

        filename = "%s/%s%.1f.png" % (self.output_folder, plane_string, depth)
        plt.gcf().set_size_inches(8, 6)
        plt.savefig(filename)

        return plot


    def get_plot_ignore_edge(self, plane_index, depth, npx, title=" ", is_normalized=False):
        plane = self.getIntensityPlane(plane_index, depth)

        xyz = [0, 1, 2]
        xyz.remove(plane_index)

        # getting number of bins in each direction
        divis = numpy.array(self.num_bins)
        # creating indices to remove edges from heat map
        #npx = 1  #  size of pixel edge to remove
        pl_ig_min1 = int( npx );  pl_ig_max1 = divis[xyz[0]] - npx
        pl_ig_min2 = int( npx );  pl_ig_max2 = divis[xyz[1]] - npx
        # truncating plane
        plane = plane[pl_ig_min1:pl_ig_max1, pl_ig_min2:pl_ig_max2]

        if is_normalized:
            plane /= numpy.max(plane.flatten())

        # truncating axes
        ax_ig_min_x = int( npx );  ax_ig_max_x = len(self.edges[0]) - npx
        ax_ig_min_y = int( npx );  ax_ig_max_y = len(self.edges[1]) - npx
        ax_ig_min_z = int( npx );  ax_ig_max_z = len(self.edges[2]) - npx
        #axes = [['x', self.edges[0]], ['y', self.edges[1]], ['z', self.edges[2]]]
        axes = [['x', self.edges[0][ax_ig_min_x:ax_ig_max_x]], ['y', self.edges[1][ax_ig_min_y:ax_ig_max_y]], ['z', self.edges[2][ax_ig_min_z:ax_ig_max_z]]]

        plane_string = axes[plane_index][0]
        del axes[plane_index]

        plt.clf()
        title = "%s = %.1f mm" % (plane_string, depth)
        plt.title(title, fontsize=22)
        plt.grid(True)
        plt.xlabel("%s [mm]" % axes[0][0], fontsize=20)
        plt.ylabel("%s [mm]" % axes[1][0], fontsize=20)

        extent = [min(axes[0][1]), max(axes[0][1]), min(axes[1][1]), max(axes[1][1])]
        plt.axis(( min(axes[0][1]), max(axes[0][1]), min(axes[1][1]), max(axes[1][1]) ))

        plot = plt.imshow(plane[::-1,:], extent=extent)

        cb = plt.colorbar(plot)
        # cb.set_ticks(ticks = [0.0,0.25, 0.5, 0.75, 1.0, 1.25])

        filename = "%s/%s%.1f.png" % (self.output_folder, plane_string, depth)
        plt.gcf().set_size_inches(8, 6)
        plt.savefig(filename)

        return plot


    def get_coordinates(self, i, j , k):

        return self.bin_centers[0][i], self.bin_centers[1][j], self.bin_centers[2][k]


    def set_maximum(self, centers, values):

        i, j, k = numpy.unravel_index(values.argmax(), values.shape)
        point = (centers[0][i], centers[1][j], centers[2][k])

        self.maximum = point

        return point


    def get_maximum(self):
        try:
            return self.maximum
        except AttributeError:
            return self.set_maximum(self.bin_centers, self.values)


    def __init__(self, input_file_path, output_folder="."):
        self.input_file = input_file_path
        self.output_folder = output_folder
        self.bins, self.edges, self.values = self.read_dose_file(self.input_file)

        self.bin_centers = self.get_bin_centers(self.edges)

        self.f3D = scipy.interpolate.RegularGridInterpolator(self.bin_centers, self.values)



if __name__ == "__main__":

    def print_usage():
        print (USAGE % (sys.argv[0], sys.argv[0]))
        sys.exit(100)

    if len(sys.argv) == 4:
        var_map = {'x':0, 'y':1, 'z':2}

        data_file_path = sys.argv[1]
        try:
            dimension = var_map[sys.argv[2].lower()]
        except KeyError:
            print_usage()

        depth = float(sys.argv[3])

    else:
        print_usage()

    my3DObj = Intensity3D(data_file_path)
    myplot = my3DObj.get_plot(dimension, depth)

    plt.show()
