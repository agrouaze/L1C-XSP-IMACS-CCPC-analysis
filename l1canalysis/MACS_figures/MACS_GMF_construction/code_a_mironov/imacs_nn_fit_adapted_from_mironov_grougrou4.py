"""
A. Grouazel
13 May 2025
this 4th version aims at integrating connexion to MLFlow and to perform training on the difference of IMACS-A(u10)*cos(phi)

"""
import pdb,glob
import logging
import time
import shutil
import numpy as np
import mlflow
from pysr import PySRRegressor
import pandas as pd
import sympy,os
from sklearn.metrics import r2_score
import argparse
import mlflow.pyfunc
from sympy import symbols, lambdify, Function
from mlflow.data import from_pandas
from datetime import datetime
import imacs_prediction_perf_figures
#tmpnpz = np.load('/home1/datahome/agrouaze/sources/git/L1C-XSP_IMACS-analysis/l1canalysis/MACS_figures/MACS_GMF_construction/code_a_mironov/imacs_interp3d_SmoothBivariateSpline.npz')
#print("tmpnpz",tmpnpz)

import pickle

# I add this class to be able to register the model trained into mlflow
class PySRModelWrapper(mlflow.pyfunc.PythonModel):
    def load_context(self, context):
        self.model = reuse_previous_model(context.artifacts["model_path"])
        # import joblib
        # self.model = joblib.load(context.artifacts["model_path"])

    def predict(self, context, model_input):
        return self.model.predict(model_input)
INPUT_DATASET_SIZE = 200000
LOWEST_WINDSPEED=3
HIGHEST_WINDSPEED=19
LOWEST_AZI=0
HIGHEST_AZI=359
LOWEST_INCIDENCE=32
HIGHEST_INCIDENCE=45
LAMBDA_MAX=100
POLARIZATION='VV'
UNIT_SAR = 'S1A'
BURSTGRP = 'intraburst'
dirpkl='/home1/datahome/agrouaze/sources/git/L1C-XSP_IMACS-analysis/l1canalysis/MACS_figures/MACS_GMF_construction/code_a_mironov/'
#with open("/home1/datahome/agrouaze/sources/git/L1C-XSP_IMACS-analysis/l1canalysis/MACS_figures/MACS_GMF_construction/code_a_mironov/imacs_interp.pkl", "rb") as f:
#with open("/home1/datahome/agrouaze/sources/git/L1C-XSP_IMACS-analysis/l1canalysis/MACS_figures/MACS_GMF_construction/code_a_mironov/imacs_interp_simplified.pkl", "rb") as f:
#with open('/home1/datahome/agrouaze/sources/git/L1C-XSP_IMACS-analysis/l1canalysis/MACS_figures/MACS_GMF_construction/code_a_mironov/imacs_interp_simplified_gaussian_20250227.pkl','rb') as f:
#with open(os.path.join(dirpkl,'imacs_interp_simplified_gaussian_macs_Im_lambda_max=100.0_limited_windrange_20250311.pkl'),'rb') as f: # with a sigma in gaussian filter that is constant=120
interpolator_file='imacs_interp_simplified_gaussian_macs_Im_lambda_max=%i.0_limited_windrange_20250321.pkl'%LAMBDA_MAX
# interpolator_file='imacs_residual_interp_simplified_gaussian_macs_Im_lambda_max=100.0_limited_windrange_20250514.pkl'
with open(os.path.join(dirpkl,interpolator_file),'rb') as f: # I used a variable sigma in gaussian filter.
    imacs_interp = pickle.load(f)
def generate_random_pts(nb_pts_sought=INPUT_DATASET_SIZE):
    #%% generate a random dataset from imacs_interp for 8000 points
    # rnd_winds = np.random.uniform(6, 15, 200000) # used for training 10 to 14
    rnd_winds = np.random.uniform(LOWEST_WINDSPEED, HIGHEST_WINDSPEED, nb_pts_sought) # for n°15 I make the hypothesis my interpolator is not that bad on a larger wind range
    rnd_azis = np.random.uniform(LOWEST_AZI, HIGHEST_AZI, nb_pts_sought)
    rnd_incs = np.random.uniform(LOWEST_INCIDENCE, HIGHEST_INCIDENCE, nb_pts_sought)
    rnd_macs = imacs_interp((rnd_winds, rnd_azis, rnd_incs))
#    rnd_residual_macs = imacs_residual_interp((rnd_winds, rnd_azis, rnd_incs))
    npz = {}
    npz['winds'] = rnd_winds
    npz['azis'] = rnd_azis
    npz['incs'] = rnd_incs
    npz['macs'] = rnd_macs
 #   npz['residual'] = rnd_residual_macs
    return npz

def createmodel(checkpoint_file=None):
    """

    method to create a new PySRRegressor but possibily reusing checkpoint from another model

    :param checkpoint_file:
    :return:
    """
    model = PySRRegressor(
        model_selection="best",  # Result is mix of simplicity+accuracy
        # niterations=4000,
        # niterations=20000,
        niterations=NB_NUMBER_ITERATIONS,
        # niterations=10,
        # temp_equation_file="equations_and_loss_IMACS.csv",  # Specify the file path
        # maxsize=35, #original
        maxsize=40,
        output_directory=new_working_dir,
        # checkpoint_file=checkpoint_file,
        populations=28,
        parsimony=1e-9,
        # weight_optimize=0.01,
        # denoise=True,
        # complexity_of_variables= 35,
        # complexity_of_operators = 20,
        # complexity_of_constants = 15,
        # precision=64,
        binary_operators=['+', '*', '/', '-', '^'],
        unary_operators=[
            "cos",
            "cos2(x) = cos(x)^2",  # Custom operator in Julia syntax
            "cos3(x) = cos(x)^3",  # Custom operator in Julia syntax
            "exp",
            "sin",
            "inv(x) = 1/x",
            "log",
            # "sm(x) = 1/(1+exp(x))",
            # "sinh",
            # "cosh",
            # "tanh",
            "tan",
            # "sec"

            # ^ Custom operator (julia syntax)
        ],
        early_stop_condition=(
            "stop_if(loss, complexity) = loss < 1e-6 && complexity < 10"
            # Stop early if we find a good and simple equation
        ),
        # constraints={"pow": (-9, 3), "mult": (-1, 20), "exp": (-1, 3), "^": (-1, 3), },
        constraints={'^': (-1, 3)},
        nested_constraints={"sin": {"sin": 0, "cos": 1}, "cos": {"sin": 1, "cos": 0}, "exp": {"exp": 1},
                            "log": {"log": 0}, "tan": {"tan": 0}, "inv": {"inv": 0}
                            },
        extra_sympy_mappings=ESM,
        # ^ Define operator for SymPy as well
        # loss="loss(x, y) = abs(x * y)"
        elementwise_loss="LogitDistLoss()"
        # ^ Custom loss function (julia syntax)
    )
    model.checkpoint_file = checkpoint_file
    #model.modeldir=new_working_dir
    return model
#OUTDIR = '/home/antoine/Documents/sources/notebooks/imacs_gmf/'
OUTDIR = '/raid/localscratch/agrouaze/imacs_gmf/iw/'
TARGET_FORM = 'residual_IMACS' # or 'IMACS'
def get_vectors():
    #%% load precalculated imac dataset from numpy file

    npz_vars = generate_random_pts()

    winds, azis, incs, imacs = npz_vars["winds"], npz_vars["azis"], npz_vars["incs"], npz_vars["macs"]
    #if TARGET_FORM=='IMACS':
    y_train = imacs*100. # for direct prediction of IMACS prediction
    #else:
    #    y_train = resi*100.
    X_train = np.array([np.deg2rad(incs), winds, np.deg2rad(azis)]).T

    #from sklearn.model_selection import train_test_split
    # here I wrote _test but in fact I use thhis part of the daatset as validation (to check that the model generalize rather than memorize)
    #X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.4, random_state=0)
    npz_vars = generate_random_pts()
    winds_val, azis_val, incs_val, imacs_val = npz_vars["winds"], npz_vars["azis"], npz_vars["incs"], npz_vars["macs"]
    #if TARGET_FORM=='IMACS':
    y_val = imacs_val*100 # for direct prediction of IMACS prediction
    #else:
    #    y_val = resi_val*100.
    X_val = np.array([np.deg2rad(incs_val), winds_val, np.deg2rad(azis_val)]).T

    weights = np.ones_like(y_train)
    weights[(azis > 0) & (azis < 50)] = 3
    weights[(azis > 310)] = 3
    weights[(azis > 150) & (azis < 220)] = 5
    return X_train,X_val,y_train,y_val,weights



ESM = {"inv": lambda x: 1 / x,
                              # "cos2": lambda x: np.cos(x)**2,
                              "cos2": lambda x: sympy.cos(x) ** 2,  # found in the doc
                              "cos3": lambda x: sympy.cos(x) ** 3,
                              # "sm": lambda x: 1 / (1+np.exp(x))
                              # "sec": lambda x: 1 / np.cos(x+np.pi)

                              }
def save_model_to_pkl(model):
    """
    method to save the pysr model into a pkl file

    :param model: pysr object
    :return:
    """
    os.makedirs(new_working_dir,exist_ok=True)
    fout = os.path.join(new_working_dir, 'pysr_fit_run_%s.pkl' % datetime.today().strftime('%Y%b%d_%H%M'))
    fid = open(fout, 'wb')
    pickle.dump(model, fid)
    fid.close()
    print('output model', model)
    logging.info('output model: %s',fout)
    # a first way to log the model in MLFLOW but column Models wont be filled
    mlflow.log_artifact(local_path=fout,artifact_path='models')
    # second way more official to register_model in MLFLOW
    # Define the artifacts and conda env
    #real_path = model.modeldir
    # real_path = model.model_["output_directory"]
    real_path = model._model._output_directory
    logging.info('real path of the model : %s',real_path)
    artifacts = {"model_path": real_path}
    npz_vars = generate_random_pts(nb_pts_sought=5)
    winds, azis, incs, imacs = npz_vars["winds"], npz_vars["azis"], npz_vars["incs"], npz_vars["macs"]
    X_sample = np.array([np.deg2rad(incs), winds, np.deg2rad(azis)]).T
    # Log and register the model
    mlflow.pyfunc.log_model(
        artifact_path="model",
        python_model=PySRModelWrapper(),
        artifacts=artifacts,
        input_example=X_sample,
        registered_model_name="PySR-%s-Model"%TARGET_FORM
    )
    return fout

def reuse_previous_model(checkpoint_file_halloffame,force_new_regressor=False):
    """

    :param checkpoint_file_halloffame: str checkpoint.pkl file path
    :return:
    """
    # line below comes from chatgpt
    #model = PySRRegressor(hall_of_fame=checkpoint_file_halloffame)

    # %% fit the previous model (version from A .Mironov may be deprecated?)
    # # model = PySRRegressor.from_file("hall_of_fame_2024-02-02_162644.080.pkl")
    #model = PySRRegressor.from_file(checkpoint_file_halloffame)
    #didi = '/home1/datahome/agrouaze/sources/git/L1C-XSP_IMACS-analysis/l1canalysis/MACS_figures/MACS_GMF_construction/code_a_mironov' 
    # logging.info('load halloffame from %s',OUTDIR)
    #lsthalloffames = sorted(glob.glob('/raid/localscratch/agrouaze/imacs_gmf/iw/202*_*_*/hall_of_fame.csv'))
    logging.info('lsthalloffames : %s',os.path.dirname(checkpoint_file_halloffame))
    if force_new_regressor is True:
        model  = createmodel(checkpoint_file=checkpoint_file_halloffame)
    else:
        model = PySRRegressor.from_file(run_directory=os.path.dirname(checkpoint_file_halloffame)) # it was working with pysr 1.5.5 (not with pysr 1.5.6)
    # model = PySRRegressor.from_file(checkpoint_file_halloffame)
    model.warm_start = True
    model.set_params(extra_sympy_mappings=ESM)
    #model.set_params(maxsize=40)
    #model.set_params(niterations=200)  # new total number of iterations
    #
    #model.warm_start = True
    #from pysr.julia_helpers import init_julia
    #
    #init_julia()
    #
    #from julia import SymbolicRegression  # Needed to load library (usually this is done by .fit())
    #from julia import Serialization
    #
    #chkptfile = os.path.join(os.path.dirname(checkpoint_file_halloffame),"checkpoint.pkl")
    #logging.info('chkptfile ; %s',chkptfile)
    #model.raw_julia_state_ = Serialization.deserialize(chkptfile)
    return model

def first_fit(metaiteration=0,outputlog_validation=None,tot_iteration=0,pklcheckpoint=None,model=None):
    """

    :param metaiteration: integer to indicate the number of training performed on the same dataset with same configuration
    :param outputlog_validation: full path of a CSV file
    :param pklcheckpoint: str path to use if you want to restart from  a previous training [default is None]
    :param model: PySRRegressor object load from a previous fit (in the same training run) [optional, default=none]
    :return:
    """
    X_train,X_val,y_train,y_val,weights = get_vectors()
    #if metaiteration!=0:
    #if model is not None or pklcheckpoint is not None:
    if model is None and pklcheckpoint is None:
        logging.info('brand new model regressor')
        model = createmodel()
    else:
        #checkpoint_file_halloffame = pklmodelsaved # here i suppose I can replace the checkpoint file by the model saved TBC
        if model is not None:
            checkpoint_file_halloffame = model.get_pkl_filename()
            logging.info('reuse the model from previous fit in the same run: %s',checkpoint_file_halloffame)
            force_new_regress = False
        else:
            checkpoint_file_halloffame = pklcheckpoint
            logging.info('pklcheckpoint : %s',pklcheckpoint)
            logging.info('reuse the model from a previous run (from a file path )')
            force_new_regress = True
        model = reuse_previous_model(checkpoint_file_halloffame,force_new_regressor=force_new_regress)
    #logging.info('set the new working dir to : %s',new_working_dir)
    #model.modeldir = new_working_dir
    model.batching = True
    model.batch_size = 300
    # model.return_state = True
    model.fit(X_train, y_train, weights=weights)
    logging.info('fit is finished.')
    #logging.info('hall of fame file : %s',model.get_hall_of_fame())
    num_iterations = len(model.equations_) # answer chatgpt to be checked, it seems weird to me
    tot_iteration += num_iterations
    # print(model)
    # Evaluate on validation set
    y_pred_val = model.predict(X_val, index=-1)
    y_pred_train = model.predict(X_train, index=-1)
    train_loss = np.mean((y_pred_train - y_train) ** 2)  # Compute validation loss, MSE
    val_loss = np.mean((y_pred_val - y_val) ** 2)  # Compute validation loss, MSE

    r2 = r2_score(y_val, y_pred_val)
    mlflow.log_metric("r2_score_validation", r2,step=metaiteration)
    mlflow.log_metric("PySR_score", model.equations_.iloc[0]['score'], step=metaiteration)
    mlflow.log_metric("PySR_score", model.equations_.iloc[0]['complexity'], step=metaiteration)

     # TODO add a log file with train loss + val loss + total iterations + score +equation most complex
    # print("Validation Loss:", val_loss)
    if outputlog_validation is not None:
        # open existing file
        if os.path.exists(outputlog_validation):
            df = pd.read_csv(outputlog_validation)
        else:
            
            df = pd.DataFrame(columns=["script", "metaiteration", "total_iterations", "val_loss", "train_loss"])
        # fid = open(outputlog_validation,'a')
        log_data = {
            'script': os.path.basename(__file__),
            'metaiteration': metaiteration,
            'total_iterations': tot_iteration,
            'val_loss': round(val_loss, 3),
            'train_loss': round(train_loss, 3)
        }
        df.loc[len(df)] = log_data
        mlflow.log_metric("val_loss", val_loss, step=metaiteration)
        mlflow.log_metric("train_loss", train_loss, step=metaiteration)

        # fid.write('MACS-pysr %s nb-train: %i, total iterations: %i val-loss: %1.3f train-loss: %1.3f'%(os.path.basename(__file__),
        #                                                                                                metaiteration,tot_iteration,val_loss,train_loss))
        df.to_csv(outputlog_validation)
        # fid.close()


        # with open(outputlog_validation, 'a') as fid:
        #     yaml.dump([log_data], fid, default_flow_style=False, allow_unicode=True)
    return tot_iteration,model

def plot_loss_training(run_name,loss_file,model):
    figpath = os.path.join(new_working_dir,'imacs_pysr_loss_evolution_%s.png'%run_name)
    import matplotlib.pyplot as plt
    # the loss given in model.equations_ will show only the latest training epochs, while I am more interested in the validation loss evolution per fit loop
    #equations_df = model.equations_
    #plt.plot(equations_df["iteration"], equations_df["loss"], marker="o", linestyle="-",label='training loss (MSE) of last fit')
    df = pd.read_csv(loss_file,names=['script','metaiteration','total_iterations','val_loss','train_loss'],header=0)
    plt.plot(df['metaiteration'],df['val_loss'],'r.-',label='val loss')
    plt.plot(df['metaiteration'],df['train_loss'],'b.-',label='train loss')
    plt.legend()
    plt.xlabel("Iteration (loop of fit)")
    plt.ylabel("Loss")
    plt.title("Loss Evolution in PySR %s"%run_name)
    plt.grid()

    plt.savefig(figpath)
    logging.info('figpath: %s',figpath)
    mlflow.log_artifact(local_path=figpath, artifact_path="figures")

    # plt.show()
def copy_pysr_working_dir(src_dir, dst_dir, overwrite=False):
    """
    Copy a full PySR working directory to a new location.

    Parameters:
    - src_dir (str): Path to the original PySR working directory.
    - dst_dir (str): Path to the new destination directory.
    - overwrite (bool): If True, overwrite dst_dir if it exists.
    """
    if not os.path.isdir(src_dir):
        raise ValueError(f"Source directory '{src_dir}' does not exist.")

    if os.path.exists(dst_dir):
        if overwrite:
            shutil.rmtree(dst_dir)
            logging.info(f"Overwriting existing directory: {dst_dir}")
        else:
            raise FileExistsError(f"Destination '{dst_dir}' already exists. Use overwrite=True to replace.")

    shutil.copytree(src_dir, dst_dir)
    logging.info(f"Copied PySR working directory from '{src_dir}' to '{dst_dir}'.")



def write_imacsAF_py(eq_str: str, filename="temporary_analytical_formula.py"):
    with open(filename, "w") as f:
        f.write("import numpy as np\n\n")
        f.write("def imacsAF(phi, alpha, u):\n")
        f.write("    \"\"\"\n")
        f.write("    :param phi: (np.ndarray) azimuth wind direction (relative to antenna), clockwise degrees x2 [radians]\n")
        f.write("    :param alpha: (np.ndarray) incidence angles x0 [degrees]\n")
        f.write("    :param u: (np.ndarray) wind speed x1\n")
        f.write("    :return: imacs (np.ndarray): GMF-like prediction (normalized)\n")
        f.write("    \"\"\"\n")
        f.write("    alpha = np.deg2rad(alpha)\n")
        f.write("    fit = " + eq_str + "\n")
        f.write("    imacs = fit / 100\n")
        f.write("    return imacs\n")

def import_the_new_method():
    import importlib.util
    import sys
    from pathlib import Path

    # Chemin vers ton fichier
    module_path = Path("temporary_analytical_formula.py")

    # Nom du module temporaire
    module_name = "temporary_analytical_formula"

    # Spécification d'import
    spec = importlib.util.spec_from_file_location(module_name, str(module_path))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = mod
    spec.loader.exec_module(mod)

    # Maintenant tu peux utiliser la fonction !
    imacsAF = mod.imacsAF
    return imacsAF

if __name__ =='__main__':
    root = logging.getLogger()
    if root.handlers:
        for handler in root.handlers:
            root.removeHandler(handler)
    parser = argparse.ArgumentParser(description="train-pysrr-imacs-gmf")
    parser.add_argument("--verbose", action="store_true", default=False)
    parser.add_argument("--model2start", action="store", default=None,
                        help='path of the checkpoint.pkl (the full working dir will be copied in a new working dir)')
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
    logging.info('start')
    # Set our tracking server uri for logging
    # mlflow.set_tracking_uri(uri="http://127.0.0.1:6080") # started on compute-101-23 http://134.246.184.23:6080
    mlflow.set_tracking_uri(uri='http://134.246.184.23:6080')
    mlflow.set_experiment("pysr_%s_experiments"%TARGET_FORM)

    # NB_NUMBER_META_ITERATIONS = 500
    NB_NUMBER_META_ITERATIONS = 5
    NB_NUMBER_ITERATIONS = 30
    logging.info('go for %i meta-iterations.',NB_NUMBER_META_ITERATIONS)
    logging.info('number of iteration per meta-iterations: %i',NB_NUMBER_ITERATIONS)
    # Check if a run is already active
    if mlflow.active_run():
        logging.info('there is active mlflow run')
        mlflow.end_run()  # Manually end it
    with mlflow.start_run():
        mlflow.log_param('script', os.path.basename(__file__))
        date_train = datetime.today().strftime('%Y%b%d-%H%M')
        tot_iteration = 0
        logging.info('model to start with: %s',args.model2start)
        run_name = 'MACS-pysr_training_loss_evolution_%s'%(date_train)
        loss_file = os.path.join(OUTDIR,run_name+'.txt')
        new_working_dir = os.path.join(OUTDIR,date_train)
        mlflow.log_param('new_working_dir',new_working_dir)
        mlflow.log_param('interpolator',interpolator_file)
        # I do a single get_vectors just to log it into MLFlow, but at each meta iteration there are new vectors generated
        X_train, X_val, y_train, y_val, weights = get_vectors()
        # Combine everything into one DataFrame
        train_df = pd.DataFrame(X_train, columns=[f"x{i}" for i in range(X_train.shape[1])])
        train_df["y"] = y_train
        train_df["weight"] = weights
        dataset = from_pandas(
            df=train_df,
            #source= 'pandas',# ("L1C-wave_training_dataset_B09","gaussian_fit",'3D-CUBE','RegularGridInterpolator'),
            targets= 'y',
            name='wave_training_dataset_B09-%s-dataset'%TARGET_FORM,
        # digest: Optional[str] = None,
        # predictions: Optional[str] = None,
        )
        mlflow.log_input(dataset, context="training")
        # end of the declaration of the datasets in mlflow

        mlflow.log_param('target',TARGET_FORM)
        if args.model2start is not None:
            copy_pysr_working_dir(os.path.dirname(args.model2start), new_working_dir)
            copied_checkpoint_file_to_start = os.path.join(new_working_dir,os.path.basename(args.model2start))
        else:
            copied_checkpoint_file_to_start = None
        t0 = time.time()
        for metaiter in range(NB_NUMBER_META_ITERATIONS):
            logging.info('\n======== run %i =======\n\n',metaiter)
            if metaiter==0:
                mlflow.log_param('model2startWith',args.model2start)
                tot_iteration,model = first_fit(metaiteration=metaiter, outputlog_validation=loss_file,
                                                       tot_iteration=tot_iteration,pklcheckpoint=copied_checkpoint_file_to_start,model=None)
            else:
                tot_iteration,model = first_fit(metaiteration=metaiter, outputlog_validation=loss_file,
                                                        tot_iteration=tot_iteration, pklcheckpoint=None,model=model)
        pklout_model = save_model_to_pkl(model)
        elapsed_time_train = time.time()-t0
        mlflow.log_params({
            "total_number_of_iterations": tot_iteration,
            "added_number_of_meta_iterations": NB_NUMBER_META_ITERATIONS,
            "number_of_iteration_per_loop": NB_NUMBER_ITERATIONS,
            "binary_operators":  model.binary_operators,
            "unary_operators": model.unary_operators,
            "model_selection": model.model_selection,
            "elapsed_time_train":elapsed_time_train,
            'dataset_size':INPUT_DATASET_SIZE,
            'highest_wind_speed':HIGHEST_WINDSPEED,
            'lowest_wind_speed':LOWEST_WINDSPEED,
            'highest_incidence': HIGHEST_INCIDENCE,
            'lowest_incidence': LOWEST_INCIDENCE,
            'highest_azimuth': HIGHEST_AZI,
            'lowest_azimuth': LOWEST_AZI,
        })
        best_equation = str(model.sympy())  # Or model.get_best()
        mlflow.log_param("best_equation", best_equation)
        logging.info('end of fit loops')
        equations_df = model.equations_
        equations_df.to_csv("equations.csv", index=False)
        mlflow.log_artifact("equations.csv",artifact_path='equations')
        # generate a temporary .py file with the function containing the analytical equation

        # 1. Extraire l'expression SymPy
        sym_expr = model.get_best(as_string=False)

        # 2. Définir les symboles avec noms explicites
        phi, alpha, u = symbols("phi alpha u")

        # 3. Création d'un mapping entre anciens noms et nouveaux
        subs_dict = {
            sympy.Symbol("x0"): alpha,
            sympy.Symbol("x1"): u,
            sympy.Symbol("x2"): phi,
        }

        # 4. Remplacement dans l'expression
        sym_expr_named = sym_expr.subs(subs_dict)

        # 5. Génération du code Python lisible
        from sympy.printing.pycode import pycode

        equation_str_numpy = pycode(sym_expr_named)
        write_imacsAF_py(equation_str_numpy)
        imacsformula = import_the_new_method()
        # plots part
        plot_loss_training(run_name,loss_file,model=pklout_model)
        fout = os.path.join(new_working_dir, 'imacs_vs_azimuth_per_windspeed_%s.png' % run_name)
        imacs_prediction_perf_figures.imacs_vs_azimuth_per_windspeed(imacs_interp,imacsformula,fout=fout)
        mlflow.log_artifact(local_path=fout, artifact_path="figures")

        fout = os.path.join(new_working_dir, 'imacs_vs_incidence_%s.png' % run_name)
        imacs_prediction_perf_figures.imacs_vs_incidence(imacs_interp,imacsformula,fout=fout)
        mlflow.log_artifact(local_path=fout, artifact_path="figures")

        fout = os.path.join(new_working_dir, 'imacs_vs_inc_12subplots_%s.png' % run_name)
        imacs_prediction_perf_figures.imacs_vs_inc_12subplots(imacs_interp, my_pola=POLARIZATION,
                                                              lambda_val=LAMBDA_MAX, burst_family=BURSTGRP, fout=fout)
        mlflow.log_artifact(local_path=fout, artifact_path="figures")

        fout = os.path.join(new_working_dir, 'imacs_vs_windspeed_fourches_%s.png' % run_name)
        imacs_prediction_perf_figures.imacs_vs_windspeed_fourches(imacsformula, my_pola=POLARIZATION,
                                                    lambda_val=LAMBDA_MAX, burst_family=BURSTGRP, sar_unit=UNIT_SAR, fout=fout)
        mlflow.log_artifact(local_path=fout, artifact_path="figures")

        fout = os.path.join(new_working_dir, 'imacs_vs_azimuth_4subplots_%s.png' % run_name)
        analyticalform_imacs_mean,interpolator_imacs_mean,obs_imacs_mean = imacs_prediction_perf_figures.imacs_vs_azimuth_4subplots(imacs_interp,
                                                        imacsformula, my_pola=POLARIZATION, lambda_val=LAMBDA_MAX,
                                                                            burst_family=BURSTGRP, sar_unit=UNIT_SAR, fout=fout)
        mlflow.log_artifact(local_path=fout, artifact_path="figures")

        imacs_prediction_perf_figures.get_table_residuals_obs_vs_interpoaltor(obs_imacs_mean, interpolator_imacs_mean)

        imacs_prediction_perf_figures.get_table_residuals_obs_vs_analyticalform(obs_imacs_mean, analyticalform_imacs_mean)
    mlflow.end_run()  # optional, but safe
