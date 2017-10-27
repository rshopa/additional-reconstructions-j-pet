# import os
import sys
import time
import numpy
import pandas
import tomopy


class ListModeReconstruction:
    """Produces object of 3D reconstruction from LOR list mode"""

    def __init__(self, geometry=None):
        """
        :param geometry: dictionary of parameters:
            N - No of detectors,
            R - PET radius, cm
            length - axial length of PET, cm
            slices - No of slices in Z direction
            C_0 - speed of light
            Delta_s, Delta_fi, Delta_z - quantization by displacement, angle
                                         and axial coordinate, respectively
            centered - if True, the Z coordinate of each LOR is estimated from its centre,
                       otherwise hit time differences are considered
        """
        if geometry is None:
            self.geometry = {'N': 384,
                             'R': 43.73,
                             'length': 50.,
                             'slices': 200,
                             'C_0': 3e8,
                             'centered': False}
        else:
            self.geometry = geometry
        self.geometry['Delta_s'] = self.geometry['R'] * numpy.pi / self.geometry['N']
        self.geometry['Delta_fi'] = numpy.pi / (self.geometry['N'] - 1)  # 0 to 199, i.e. 0 to N-1
        self.geometry['Delta_z'] = self.geometry['length'] / (self.geometry['slices'] - 1)  # 200 points from 0 to 199

        self.lor_data = None
        self.lor_intensities = None
        self.radiograph = None
        self.reconstructed = {}

        self.__data_state = {        # flags for validation in functions
            'loaded': False,
            'mapped': False,
            'radiograph': False
        }

    def vectors(self):
        """
        Calculates axis vectors corresponding to reconstructed images,
        no use as a separate object attribute is needed.
        :return: dict of 1-dimensional arrays
        """
        # X and Y axis
        xy_points = int(round(self.geometry['R'] / self.geometry['Delta_s'])) * 2 + 1
        xy = numpy.subtract(numpy.multiply(self.geometry['Delta_s'], range(xy_points)),
                            self.geometry['Delta_s'] * (xy_points - 1) / 2)
        # Z axis
        z = numpy.subtract(numpy.multiply(self.geometry['Delta_z'], range(self.geometry['slices'])),
                           self.geometry['Delta_z'] * (self.geometry['slices'] - 1) / 2)
        return {'x': xy, 'y': xy, 'z': z}

    def read_data(self, filename):
        """
        impors data from file in list mode
        :param filename: path to file
        :return: none
        """
        sys.stdout.write("Loading data... ")
        self.lor_data = pandas.read_table(filename, header=None)
        sys.stdout.write("Done!\nData has been successfully imported.\n")
        sys.stdout.flush()
        self.__data_state['loaded'] = True  # state change

    def to_indices_with_intensities(self, dat=None):
        """
        Creates dataframe with radiographs' (or sinograms?)
        :param dat: dataframe of lors
        :return: dataframe with indices and intensities
        """
        if not self.__data_state['loaded']:
            print "ERROR: No data founded. Please, import your data first."
            return "Aborted."
        if dat is None:
            dat = self.lor_data  # default value is the whole set
        lor_norm = pandas.DataFrame({'fi_index': [],
                                     'z_index': [],
                                     's_index': []}, dtype=int)  # empty dataframe
        sys.stdout.write('Mapping data:')  # setup progress bar
        # iterate
        for index, row in dat.iterrows():
            new_lor = self.__map_to_indices(self.__rebinned_lor(row))  # read one row and map to indices
            lor_norm = lor_norm.append(pandas.DataFrame({'fi_index': [new_lor['fi_index']],
                                                         'z_index': [new_lor['z_index']],
                                                         's_index': [new_lor['s_index']]}))
            step = round(float(dat.shape[0]) / 100)  # 100 dots of progress
            if index % step == 0:
                sys.stdout.write('.')  # progress dots
                sys.stdout.flush()

        sys.stdout.write(" Done!\n")
        time.sleep(0.4)
        # rescale to intensities (create additional column named 0 and rename to 'intensity')
        sys.stdout.write("Rescaling:......")
        sys.stdout.flush()
        lor_norm = lor_norm.groupby(['fi_index',
                                     's_index',
                                     'z_index']).size().reset_index().rename(columns={0: 'intensity'})
        sys.stdout.write(" Done!\n")
        sys.stdout.flush()
        self.__data_state['mapped'] = True  # change 'mapped' state into True
        return lor_norm

    def create_radiograph(self):
        """
        Creates 3-dimensional array of intensities
        :return: none
        """
        if not self.__data_state['mapped']:
            print "WARNING: Data has not been mapped yet.\n"
            if self.__query_yes_no("Do you want to map it now?"):
                self.lor_intensities = self.to_indices_with_intensities()  # map initial data
                if type(self.lor_intensities) == str:
                    self.lor_intensities = None
                    return "Aborted!"

            else:
                return "Aborted."
        # initialize empty array
        sys.stdout.write("Creating radiograph...")
        sys.stdout.flush()
        self.radiograph = numpy.zeros((self.geometry['N'],
                                       self.geometry['slices'],
                                       int(2 * self.geometry['R'] / self.geometry['Delta_s']) + 1), dtype=int)
        # assign intensities
        for index, row in self.lor_intensities.iterrows():
            self.radiograph[row['fi_index'], row['z_index'], row['s_index']] = row['intensity']

        sys.stdout.write("Done!\n")
        sys.stdout.flush()
        self.__data_state['radiograph'] = True  # self state change

    def add_rectonstruction(self, filter=None, projections=None):
        """
        Appends reconstructed image to list self.reconstructed (3-dimensional array)
        :param filter: 'none', 'fbp' (ASTRA FBP) are available
        :param projections: No of angles
        :return: info message
        """
        if projections is None:
            projections = self.geometry['N']   # default
        __angles = tomopy.angles(projections)  # internal array parameter of tomopy
        if not self.__data_state['radiograph']:
            print "ERROR: No radiograph founded!"
            return "Aborted."

        if filter is None or filter.lower() == 'none':
            sys.stdout.write("Reconstructing without filtering... ")
            sys.stdout.flush()
            # NOTE! Since no tomopy fbp implementation has been scripted yet,
            # it is used for the reconstruction without filtering
            self.reconstructed['none'] = tomopy.recon(self.radiograph, __angles,
                                                      algorithm='fbp', filter_name='none')
        elif filter == 'fbp':
            sys.stdout.write("Reconstructing with ASTRA FBP... ")
            sys.stdout.flush()
            # ASTRA algorithm is used for the FBP reconstruction
            __options = {'proj_type': 'linear', 'method': 'FBP'}
            self.reconstructed['fbp'] = tomopy.recon(self.radiograph, __angles,
                                                     algorithm=tomopy.astra, options=__options)
        else:
            print "ERROR: Unsupported filter type"
            return "Aborted."

        sys.stdout.write("Done!\n")
        sys.stdout.flush()
        return "Successfully inserted into dictionary"

    def export_slice(self, file_name, axis, index, filt='none', vect=True, sep=' '):
        """
        Exports slice into file
        :param file_name: where to export
        :param axis: name of axis, i.e. 'x', 'y, 'z'
        :param index: index of slice
        :param filt: key to dictionary of reconstructed images (filter name)
        :param vect: whether to export vectors of axis values
        :param sep: delimiter
        :return: saves 2D slice into file
        """
        axis = axis.lower()
        if axis.lower() not in 'xyz':
            print "Inappropriate index selected."
            return "Aborted."
        try:
            img = self.reconstructed[filt]
        except KeyError:
            print "ERROR: No such image"
            return "Aborted."
        # select slice (??? how to avoid multiple ifs ???)
        if axis == 'x':
            cut = img[:, :, index]
        elif axis == 'y':
            cut = img[:, index, :]
        else:
            cut = img[index, :, :]
        # save to file
        vector_keys = [k for k in 'xyz' if k != axis]   # exclude axis of slice cut
        prefix = "".join(sorted(vector_keys))           # will be prepended to file name
        numpy.savetxt(prefix + "_" + filt + "_" + file_name, cut, delimiter=sep)
        # save separate files for axis of the slice if vect selected
        if vect:
            for k in vector_keys:
                # vector_file_name = k+"_axis_"+filt+"_"+file_name  # maybe no filter should be specified
                # vector_file_name = k + "_axis_" + file_name
                vector_file_name = k + "_axis"
                numpy.savetxt(vector_file_name, self.vectors()[k], delimiter=sep)
        print "Successfully saved to file " + prefix + "_" + filt + "_" + file_name
        return "Done!"

    # PRIVATE METHODS ########################################

    def __centre_lor(self, dat):
        """
        calculates centre of LOR
        :type dat: pandas dataframe
        :param dat: single row of listmode dataframe
        :return: dictionary of coordinates
        """
        x_c = (dat.iloc[0] + dat.iloc[4]) / 2.
        y_c = (dat.iloc[1] + dat.iloc[5]) / 2.
        z_c = (dat.iloc[2] + dat.iloc[6]) / 2.
        return dict(x=x_c, y=y_c, z=z_c)

    def __z_corrected(self, dat):
        """
        # adjusts Z coordinate according to difference of hit times
        :type dat: pandas dataframe
        :param dat: single row of listmode dataframe
        :return: adjusted coordinate Z
        """
        # multiply by 100 to obtain centimeters
        delta_r = -100 * self.geometry['C_0'] * (
        dat.iloc[3] - dat.iloc[7]) * 1e-12  # minus: inverse direction to delta_t
        # spherical coordinates:
        theta = numpy.arccos((dat.iloc[2] - dat.iloc[6]) /
                             numpy.sqrt((dat.iloc[0] - dat.iloc[4]) ** 2. +
                                        (dat.iloc[1] - dat.iloc[5]) ** 2. + (dat.iloc[2] - dat.iloc[6]) ** 2.))
        return self.__centre_lor(dat)['z'] + (delta_r * numpy.cos(theta)) / 2.

    def __rebinned_lor(self, row):
        """
        rebin LOR to TomoPy like list - NOT ROUNDED!
        :type row: pandas dataframe
        :param row: single row of listmode data frame
        :return: dictionary of 3d sinogram {phi,z,s}
        """
        centres = self.__centre_lor(row)
        if self.geometry['centered']:
            z = centres['z']  # no change
        else:
            z = self.__z_corrected(row)  # adjust

        # negative displacement if y is below zero
        displacement = numpy.sqrt(centres['x'] ** 2 + centres['y'] ** 2) * \
                       (numpy.sign(centres['y']) + (numpy.sign(centres['y']) == 0) * numpy.sign(centres['x']))
        try:
            fi = numpy.arctan(-float(row[4] - row[0]) / float(row[5] - row[1]))  # a trick to allow ZeroDivisionError
        except ZeroDivisionError:
            fi = numpy.pi / 2.

        # if arctan(fi)<0 - add +pi
        angle = fi + numpy.pi * (fi < 0)
        return dict(angle=angle, z=z, displacement=displacement)  # dictionaries are 'faster' than lists

    def __map_to_indices(self, lor, maximums=None):
        """
        returns indices for TomoPy 3d reconstruction data
        :param lor: rebinned lor
        :param maximums: maximum values of each parameter
        :return: dictionary of 3D array indices
        """
        if maximums is None:
            maximums = {'angle': numpy.pi,
                        'z': self.geometry['length'] / 2,
                        'displacement': self.geometry['R']}

        angle_index = int(round(lor['angle'] / self.geometry['Delta_fi']))
        z_index = int(round((lor['z'] + maximums['z']) / self.geometry['Delta_z']))
        displacement_index = int(round((lor['displacement'] +
                                        maximums['displacement']) / self.geometry['Delta_s']))
        return dict(fi_index=angle_index, z_index=z_index, s_index=displacement_index)

    def __query_yes_no(self, question, default="yes"):
        """Ask a yes/no question via raw_input() and return their answer.

        "question" is a string that is presented to the user.
        "default" is the presumed answer if the user just hits <Enter>.
            It must be "yes" (the default), "no" or None (meaning
            an answer is required of the user).

        The "answer" return value is True for "yes" or False for "no".
        """
        valid = {"yes": True, "y": True, "ye": True,
                 "no": False, "n": False}
        if default is None:
            prompt = " [y/n] "
        elif default == "yes":
            prompt = " [Y/n] "
        elif default == "no":
            prompt = " [y/N] "
        else:
            raise ValueError("invalid default answer: '%s'" % default)

        while True:
            sys.stdout.write(question + prompt)
            choice = raw_input().lower()
            if default is not None and choice == '':
                return valid[default]
            elif choice in valid:
                return valid[choice]
            else:
                sys.stdout.write("Please respond with 'yes' or 'no' "
                                 "(or 'y' or 'n').\n")

    # def __switcher(self, axis_name):
    #     """
    #     returns index of selected axis
    #     :param axis_name: 'x', 'y' or 'z'
    #     :return: index (0, 1 or 2) or error message
    #     """
    #     try:
    #         ind = 'xyz'.index(axis_name.lower())
    #     except ValueError:
    #         print "ERROR: undefined axis"
    #         return "Aborted."
    #     return ind
