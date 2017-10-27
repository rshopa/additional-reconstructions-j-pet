import sys
import ListModeReconstruction

# Do not modify below this line
# =============================
if __name__ == '__main__':
    try:
        file_name = sys.argv[1]
    except IndexError:
        sys.exit("\nUsage: python get_slices.py <input_file> [biased] [centered]")

    sys.stdout.write("-------------------------------------\n")
    sys.stdout.flush()
    lm = ListModeReconstruction.ListModeReconstruction()
    lm.read_data(file_name)  # import data

    options = map(lambda x: x.lower(), sys.argv[2:])  # all to lowercase
    # create output name (as a subset of initial filename - last few characters)
    if 'biased' in options:
        output_name = file_name[-12:]  # x10_y0_z1875'
        z_index = int(round((lm.geometry['slices'] / 2) * (1 + .75) - 1))  # 174 as for 18.75 cm
    else:
        output_name = file_name[-9:]   # x10_y0_z0'
        z_index = int(round((lm.geometry['slices'] / 2) - 1))  # 99 for 0 cm
    if 'centered' in options:
        lm.geometry['centered'] = True
        output_name += "_centered_Z"  # append to output file name

    lm.lor_intensities = lm.to_indices_with_intensities()  # perform mapping (ALL)
    # lm.lor_intensities = lm.to_indices_with_intensities(lm.lor_data.sample(1000))  # perform mapping (sample)
    lm.create_radiograph()
    # add reconstructions
    lm.add_rectonstruction('fbp')
    lm.add_rectonstruction('none')

    # export slices to files ('novector to avoid duplicating')
    if 'novectors' in options:
        vector = False
    else:
        vector = True
    lm.export_slice(output_name, 'z', z_index, 'fbp', vect=vector)
    lm.export_slice(output_name, 'z', z_index, 'none', vect=False)  # vectors have already been created
    lm.export_slice(output_name, 'y', 122, 'fbp', vect=vector)
    lm.export_slice(output_name, 'y', 122, 'none', vect=False)
