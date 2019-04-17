numColWhereRedshift = 2
numBins = 12
lastNumberBin = numBins - 1
maxRedshift = 0.30
distinctionRedshiftInBin = maxRedshift / numBins
deltaOne = 0.19
deltaTwo = 0.15
offsetAverage = 0.03
numClusters = 2712
minNumberOfGalaxiesInFile = 15

class GalaxyRecord:
    def __init__(self, str_record):
        self.record = str_record
        self.g = float(str_record.split()[numColWhereRedshift])


def read_data_from_file_with_cluster(file_name):
    fd_input = open(file_name, 'r')
    list_gal = list()
    for data_str in fd_input.readlines():
        list_gal.append(GalaxyRecord(data_str))
    return list_gal


def number_of_bin_galaxy_belong(g):
    for i in range(1, numBins+1):
        end_of_bin = i * distinctionRedshiftInBin
        if g < end_of_bin:
            return i-1
    return -1
    # raise Exception


def make_distribution(list_of_gal):
    histogram = list()
    for ie in range(numBins): ##
        histogram.append(list())

    for galaxy in list_of_gal:
        number = number_of_bin_galaxy_belong(galaxy.g)
        if number != -1:
            histogram[number].append(galaxy)
    distribution_in_bins = list()

    ie = 1
    for bin_ in histogram:
        average_red_shift_in_bin = (ie - 0.5) * distinctionRedshiftInBin
        # num of galaxies in bin is O(g^2) because of cone's square
        distribution_in_bins.append(
            round(len(bin_) / average_red_shift_in_bin ** 2))
    return histogram, distribution_in_bins


def find_ind_of_max(distribution_in_bins):
    ind = 0
    i_max_0 = i_max_1 = 0
    max_0 = max_1 = 0
    for bin_ in distribution_in_bins:
        if bin_ > max_0:
            max_1 = max_0
            max_0 = bin_
            i_max_1 = i_max_0
            i_max_0 = ind
        elif bin_ > max_1:
            max_1 = bin_
            i_max_1 = ind
        ind += 1
    return i_max_0, i_max_1


def choose_bin(distribution_in_bins):
    code = 0b00000

    ind_max0, ind_max1 = find_ind_of_max(distribution_in_bins)
    bins_are_near = 1 if abs(ind_max0 - ind_max1) == 1 else 0
    all_values = 0
    for bin_ in distribution_in_bins:
        all_values += bin_
    d_m0 = distribution_in_bins[ind_max0] / all_values
    d_m1 = distribution_in_bins[ind_max1] / all_values

    print (round(d_m0, 3), round(d_m1, 3))

    if d_m0 >= deltaTwo and d_m1 < deltaTwo:
        if d_m0 >= deltaOne:
            code = 0b00000
        else:
            code = 0b00001
        return ind_max0, code

    elif d_m0 >= deltaTwo and d_m1 >= deltaTwo:
        if ind_max0 < ind_max1:
            if bins_are_near:
                code = 0b00010
                return (ind_max0, ind_max1), code
            else:
                if d_m0 >= deltaOne:
                    code = 0b00000
                else:
                    code = 0b00001
                return ind_max0, code

        elif ind_max1 < ind_max0:
            if bins_are_near:
                code = 0b00010
                return (ind_max0, ind_max1), code
            else:
                if ind_max0 == numBins - 1:
                    if d_m1 >= deltaOne:
                        code = 0b10100
                    else:
                        code = 0b00101
                    return ind_max1, code

                elif d_m1 >= deltaOne:
                    #  ind_max1 is previous and has enough galaxies
                    code = 0b10000
                    return ind_max1, code

                elif d_m1 < deltaOne:
                    if d_m0 > deltaOne:
                        code = 0b01000
                        return ind_max0, code
                    else:
                        code = 0b10001
                        return ind_max1, code

    code = 0b11111
    return ind_max0, code


def choose_galaxies_by_bin(list_of_gal, histogram, number_of_bins):
    left_bound, right_bound = detect_bounds_of_values_bins(number_of_bins)
    sum_redshift = 0
    number_of_gal_in_bounds = 0
    for i in range(left_bound, right_bound + 1):
        number_of_gal_in_bounds += len(histogram[i])
        for gal in histogram[i]:
            sum_redshift += gal.g
    average_redshift = sum_redshift / number_of_gal_in_bounds
    result_galaxies = list()
    for gal in list_of_gal:
        if ((average_redshift - offsetAverage <= gal.g) and
                (gal.g <= average_redshift + offsetAverage)):
            result_galaxies.append(gal)
    return result_galaxies, average_redshift


def detect_bounds_of_values_bins(number_of_bins):
    if type(number_of_bins) == tuple:
        less_number = min(number_of_bins)
        more_number = max(number_of_bins)
    else:
        less_number = more_number = number_of_bins

    left_bound = less_number - 1 if less_number != 0 else 0
    right_bound = more_number + 1 if more_number != lastNumberBin \
            else lastNumberBin
    return left_bound, right_bound


def record_chosen_galaxies(list_of_gal, code, average_redshift, cluster_num):
    file_name = "Abell_" + str(cluster_num) + ' ' + "less 0.3 chosen.txt"
    file_out = open(file_name, 'w')
    extra_info = str(round(average_redshift, 5)) + '\t' + code[2:].zfill(5)
    file_out.write(extra_info + '\n')
    for gal in list_of_gal:
        file_out.write(gal.record)

    file_out.close()


for i in range(0, numClusters):
    file_name = "Abell less 0.3\Abell_g_less_0.3" + str(i) + ".txt"
    list_of_gal = read_data_from_file_with_cluster(file_name)
    if len(list_of_gal) < minNumberOfGalaxiesInFile:
        continue
    histogram_all_galaxies, distribution = make_distribution(list_of_gal)
    chosen_bin, code = choose_bin(distribution)
    code = bin(code)
    list_of_chosen_galaxies, average_redshift = choose_galaxies_by_bin(
        list_of_gal, histogram_all_galaxies, chosen_bin)
    if len(list_of_chosen_galaxies) < minNumberOfGalaxiesInFile:
        continue
    record_chosen_galaxies(list_of_chosen_galaxies, code, average_redshift, i+1)
    print (len(list_of_chosen_galaxies))


