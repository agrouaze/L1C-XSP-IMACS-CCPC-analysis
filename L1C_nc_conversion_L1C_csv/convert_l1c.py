from l1c_B07_conversion import L1CConverter
import argparse

if __name__ == "__main__":
    
    # Get input parameters
    parser = argparse.ArgumentParser(description = 'Extract data of interest from L1C files and store them to pandas dataframe.')
    parser.add_argument('--input_path')
    
    args = parser.parse_args()
    input_path = args.input_path

    root_savepath = '/home1/datawork/ljessel/l1c_converted/B07'
    selected_vars = ['sigma0_filt','incidence','nesz', 'ground_heading',
                 'land_flag', 'uwnd', 'vwnd','U10', 'V10', 'macs_Im', 'macs_Re', 'normalized_variance_filt']

    # INTRABURST
    l1c_converter = L1CConverter(input_path, root_savepath, selected_vars, 'intraburst')
    l1c_converter.converter(save=True)
    
    # INTERBURST
    l1c_converter = L1CConverter(input_path, root_savepath, selected_vars, 'interburst') 
    l1c_converter.converter(save=True)
