import argparse

import ROOT

# Create an argument parser
parser = argparse.ArgumentParser(description="Run the fitting part of the analysis.")

# Add argument for the first parameter (n-max-files-per-sample)
parser.add_argument('--n-files-per-sample', type=int, required=True, help="Maximum number of files per sample.")


def get_fit_result(file_path, fit_result_name):
    """Open the ROOT file and retrieve the RooFitResult object."""
    file = ROOT.TFile(file_path)
    fit_result = file.Get(fit_result_name)
    if not fit_result:
        raise ValueError(
            f"Fit result '{fit_result_name}' not found in {file_path}"
        )
    return fit_result


def compare_fit_results(result1, result2):
    """Compare the parameter values of two RooFitResults."""
    params1 = result1.floatParsFinal()
    params2 = result2.floatParsFinal()

    # Check for the same number of parameters
    if params1.getSize() != params2.getSize():
        print(
            f"Number of parameters differ: {params1.getSize()} != {params2.getSize()}"
        )
        return

    print("Comparing parameters...")
    
    ERROR = False

    # Loop over parameters in the first result and compare with the second
    for i in range(params1.getSize()):
        par1 = params1[i]
        par2 = params2.find(
            par1.GetName()
        )  # Find corresponding parameter by name in result2

        if not par2:
            print(
                f"Parameter '{par1.GetName()}' not found in the second fit result."
            )
            ERROR = True
            continue

        # Compare values and print differences
        if abs(par1.getVal() - par2.getVal()) < 1e-6:
            print(f"Parameter '{par1.GetName()}' matches: {par1.getVal()}")
        else:
            print(
                f"Parameter '{par1.GetName()}' differs: {par1.getVal()} != {par2.getVal()}"
            )
            ERROR = True

        # Optionally compare errors too
        if abs(par1.getError() - par2.getError()) > 1e-6:
            print(
                f"Parameter '{par1.GetName()}' error differs: {par1.getError()} != {par2.getError()}"
            )
            ERROR = True
            
    if ERROR:
        print("ERROR: Comparison failed.")


args = parser.parse_args()

number_of_files = args.n_files_per_sample

# Replace these with the paths to your .root files and fit result names
file1 = "./fitResults.root"
file2 = f"./reference/fitResults/fitResults_{number_of_files}_file.root"
fit_result_name_1 = "fitResult"  # Fit result in first file
fit_result_name_2 = "fitResult"  # Fit result in second file

# Load the fit results from the two files
fit_result_1 = get_fit_result(file1, fit_result_name_1)
fit_result_2 = get_fit_result(file2, fit_result_name_2)

# Compare the fit results
compare_fit_results(fit_result_1, fit_result_2)
