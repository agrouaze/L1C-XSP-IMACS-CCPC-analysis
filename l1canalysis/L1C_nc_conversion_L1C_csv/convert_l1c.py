from l1canalysis.L1C_nc_conversion_L1C_csv.l1c_B07_conversion import L1CConverter
import argparse
import logging
if __name__ == "__main__":
    
    # Get input parameters

    parser = argparse.ArgumentParser(description = 'Extract data of interest from L1C files and store them to pandas dataframe.')
    parser.add_argument('--input_path')
    parser.add_argument("--verbose", action="store_true", default=False)
    args = parser.parse_args()

    fmt = "%(asctime)s %(levelname)s %(filename)s(%(lineno)d) %(message)s"
    if args.verbose:
        logging.basicConfig(
            level=logging.DEBUG, format=fmt, datefmt="%d/%m/%Y %H:%M:%S", force=True
        )
    else:
        logging.basicConfig(
            level=logging.INFO, format=fmt, datefmt="%d/%m/%Y %H:%M:%S", force=True
        )

    input_path = args.input_path

    # root_savepath = '/home1/datawork/ljessel/l1c_converted/B07'
    # root_savepath = '/home1/scratch/agrouaze/l1c_converted/B07extension' # used for extension dataset only
    # root_savepath = '/home1/scratch/agrouaze/l1c_converted/B07full' # used for extension dataset only
    # root_savepath = '/home1/scratch/agrouaze/l1c_converted/B07medium' # used for extension dataset only
   # root_savepath = '/home1/scratch/agrouaze/l1c_converted/B07full_with_ccpc'
    root_savepath = '/home1/scratch/agrouaze/l1c_converted/B09full_imacs_normalized'
    logging.info('root_savepath = %s',root_savepath)
    selected_vars = ['sigma0_filt','incidence','nesz', 'ground_heading',
                 'land_flag', 'uwnd', 'vwnd','U10', 'V10', 'macs_Im', 'macs_Re', 'normalized_variance_filt',
                 'CCPC_filt_Re','CCPC_filt_Im']
    # add variables from WW3
    selected_vars += ['t0m1','hs']
    # INTRABURST
    l1c_converter = L1CConverter(input_path, root_savepath, selected_vars, 'intraburst')
    l1c_converter.converter(save=True)


    # INTERBURST
    selected_vars = ['sigma0_filt','incidence','nesz', 'ground_heading',
                 'land_flag', 'uwnd', 'vwnd','U10', 'V10', 'macs_Im', 'macs_Re', 'normalized_variance_filt',
                 'CCPC_overlap_filt_Re','CCPC_overlap_filt_Im']
    l1c_converter = L1CConverter(input_path, root_savepath, selected_vars, 'interburst') 
    l1c_converter.converter(save=True)
    logging.info('end of processing : success')
